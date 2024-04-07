/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "Colvar.h"

#include "tools/Matrix.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>
#include <iostream>
using namespace std;

namespace PLMD {
namespace colvar {

class VoronoiC0 : public Colvar {
  bool pbc;
  bool serial;
  std::unique_ptr<NeighborList> nl;
  std::vector<PLMD::AtomNumber> list_a,list_b,list_c;
  std::vector<PLMD::AtomNumber> atomsToRequest;
  bool invalidateList;
  bool firsttime;
  int  lambda;
  //int groupasel;
  int nrx, num_atomsa,num_atomsb,num_atoms,num_atomso;
  double d0, d1, d2, d3, r0, sum_exp;

public:
  explicit VoronoiC0(const ActionOptions&);
  ~VoronoiC0();
// active methods:
  void calculate() override;
  void prepare() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(VoronoiC0,"VORONOIC0")

void VoronoiC0::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbor list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbor list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbor list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("compulsory","LAMBDA","1","The lambda parameter of the sum_exp function; 0 implies 1");  
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","D_1","0.0","The d_1 parameter of the switching function");
  keys.add("compulsory","D_2","0.0","The d_2 parameter of the switching function");
  keys.add("compulsory","D_3","0.0","The d_3 parameter of the switching function");
  keys.add("compulsory","NRX","0.0","The number of reactive sites");
  //keys.add("compulsory","GROUPASEL","1","Number of atoms in groupA to be used in the CV");
}

VoronoiC0::VoronoiC0(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  invalidateList(true),
  firsttime(true)
{

  parseFlag("SERIAL",serial);

  std::vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);
  
  //parse("GROUPASEL",groupasel);

  list_a = ga_lista;
  list_b = gb_lista;
  
  num_atomsa = list_a.size();
  num_atomsb = list_b.size(); 
  num_atoms = num_atomsa + num_atomsb;  

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;  
  
  parse("D_0",d0);
  parse("D_1",d1);
  parse("D_2",d2);
  parse("D_3",d3);
  parse("NRX",nrx);
  parse("LAMBDA",lambda);  
  num_atomso = num_atomsa - nrx;  

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh) {
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }

  addValueWithDerivatives(); setNotPeriodic();
  if(gb_lista.size()>0) {
    if(doneigh)  nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm,nl_cut,nl_st);
    else         nl=Tools::make_unique<NeighborList>(ga_lista,gb_lista,serial,dopair,pbc,getPbc(),comm);
  } else {
    if(doneigh)  nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm,nl_cut,nl_st);
    else         nl=Tools::make_unique<NeighborList>(ga_lista,serial,pbc,getPbc(),comm);
  }

  requestAtoms(nl->getFullAtomList()); 

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0; i<gb_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }
}

VoronoiC0::~VoronoiC0() {
// destructor required to delete forward declared class
}

void VoronoiC0::prepare() {
  if(nl->getStride()>0) {
    if(firsttime || (getStep()%nl->getStride()==0)) {
      requestAtoms(nl->getFullAtomList()); 
      invalidateList=true;
      firsttime=false;
    } else {
      //requestAtoms(nl->getReducedAtomList());
      requestAtoms(nl->getFullAtomList()); //I had to take the FULL list all the time, otherwise the different order in the local list of atoms create a lot of trouble later. Slightly more inefficient but unavoidable      
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void VoronoiC0::calculate()
{

  double totcharge=0.0;
  vector<double> sum_exp(getNumberOfAtoms());
  fill(sum_exp.begin(),sum_exp.end(),0.);

  //getNumberOfAtoms() gives the total number of atoms requested. It depends on Nlist. Its index is LOCAL, not absolute!!! Every stride it is ALL atoms requested and the order is different.
  //the link between local and absolute is getAbsoluteIndexes()[k]. It is a LOCAL list, where the k local element indicates its ABSOLUTE index !!!!!!!!!!
  //nn below is different, it's the list of couples of atoms within range

  Tensor virial;
  vector<Vector> deriv(getNumberOfAtoms());
  //Vector zeros;
  //zeros.zero();
  //fill(deriv.begin(), deriv.end(), zeros);


  if(nl->getStride()>0 && invalidateList) {
    nl->update(getPositions());
  }

  unsigned stride;
  unsigned rank;
  if(serial) {
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned nt=OpenMP::getNumThreads();
  const unsigned nn=nl->size(); //it is the list of neighbors!
  //first, atom0 in group A and all its neighbours from groupB, then another from A and so on. It is a way to put all couples in a 1D list
  if(nt*stride*10>nn) nt=1;

  //NEW INSTRUCTIONS!!!!!!
  //vector<Vector> nndistvect(nn);
  //vector<double> nndist(nn);
  vector<double> nnexp(nn);
  vector<unsigned> nni0(nn);
  vector<unsigned> nni1(nn);
  vector<double> nnexpnorm(getNumberOfAtoms()); //normalization of the colvar, local list
  vector<vector<double>> c(getNumberOfAtoms(),vector<double>(getNumberOfAtoms()));  //store exponentials in local indeces
  vector<double> charge(num_atomsa); //if I request reducedatomlist, it is a lot of trouble here because the local index of atoms in groupa changes in the update step
  vector<vector<Vector>> distAB(num_atomsa,vector<Vector>(getNumberOfAtoms())); //I store all distances between groupA and B so not to compute them again in the derivative loop
  vector<vector<double>> distABinvmod(num_atomsa,vector<double>(getNumberOfAtoms())); //the inverse modulus of distAB, it saves time to store it and read it in the double loop
  
  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(getPositions().size());  //defined locally so that each thread will modify it independently. Later all will gather
    vector<double> charget(num_atomsa); //these are defined locally to each thread so the loops don't get stopped all the time. Then their global counterpart is updated
    vector<double> nnexpnormt(getNumberOfAtoms());
    //LOOP1 compute the distances, the exponentials and the normalizations over groupb elements
    #pragma omp for
    for(unsigned int i=0; i<nn; i+=1) {

        unsigned i0=nl->getClosePair(i).first;  //local index
        unsigned i1=nl->getClosePair(i).second;

        //if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;

        if(pbc) {
          distAB[i0][i1]=pbcDistance(getPosition(i0),getPosition(i1));
        } else {
          distAB[i0][i1]=delta(getPosition(i0),getPosition(i1));
        }

        distABinvmod[i0][i1]=1.0/distAB[i0][i1].modulo();
        nnexp[i]=exp(lambda * distAB[i0][i1].modulo());
        nni0[i]=i0;  //local index of the atom in groupA
        nni1[i]=i1;
        //printf("%d %f %e \n",i,nndist[i],nnexp[i] );  
        //printf("%d %d %d \n",i,nn,getNumberOfAtoms() );
        //printf("local %d %d couple %d \n", i0,i1,i); 
        //printf("absolute index %d %d couple %d \n", getAbsoluteIndex(i0),getAbsoluteIndex(i1),i);    

        //to fill this, I need the LOCAL index in the reduced list, not i1!!
        //#pragma omp critical
        nnexpnormt[i1]+=nnexp[i]; //build normalization, only on second index
        //printf("%d    %d %d %lf   %le %le \n",i,nni0[i],nni1[i],nndist[i],nnexp[i],nnexpnorm[nni1[i]]);   
    }
    #pragma omp critical
    for(unsigned i=0; i<getNumberOfAtoms(); i++) nnexpnorm[i]+=nnexpnormt[i];
    #pragma omp barrier

    //LOOP2 compute the unshifted charge on atoms A
    #pragma omp for
    for(unsigned int i=0; i<nn; i+=1) {
      c[nni0[i]][nni1[i]]=nnexp[i]/nnexpnorm[nni1[i]];  //the couple of indeces is unique, there is no overwriting in the array
      charget[nni0[i]]+=c[nni0[i]][nni1[i]];
    }
    #pragma omp critical
    for(unsigned i=0; i<num_atomsa; i++) charge[i]+=charget[i];
    //for(unsigned i=0; i<groupasel; i++) charge[i]+=charget[i]; 
    //for(unsigned i=0; i<num_atomso; i++) charge[i]+=charget[i];  
    #pragma omp barrier

    //LOOP2.5 shift charge by d0
    //calculate the charge of H2O. -1 is negative, 0 is neutral, and 1 is positive
    #pragma omp for reduction(+:totcharge)
    for(unsigned int j=0; j<num_atomso; j+=1) {    
      charge[j]-=d0;//d1 d2 d3 is not used here
      //totcharge+=charge[j]; 
      totcharge+=pow(charge[j],2);	  
    }    

    //LOOP4 compute the derivatives on couples AB. OPTIMIZE double loop
    //further loop on the atoms of groupA, but I do it on the couples    
    //COSTLY STEP, CHECK IT  
    //vector<vector<Vector>> derivt(nn,vector<Vector>(num_atomsa));  //derivative local to each thread, same indeces
    #pragma omp for collapse(2)
    for(unsigned int i=0; i<nn; i+=1) {  //it is equivalent to the loop on n,j of Emanuele. The derivative will be evaulated on these two atoms
        //for(unsigned int k=0;k<num_atomsa;k++) {
        //for(unsigned int k=0;k<groupasel;k++) {
		for(unsigned int k=0; k<num_atomsa; k+=1) {  
        int d_in,ind0,ind1;
        double buf,fk;
        //Vector distance_kind1;

        ind0=nni0[i]; //look for the groupA index related to the couple of atoms in i
        ind1=nni1[i]; //same for the groupB index
		
		if (ind0 == k) d_in = 1;
		if (ind0 != k) d_in = 0;
		fk = lambda * c[k][ind1] * (d_in - c[ind0][ind1]);
		
		if(ind0 > num_atomso-1) {
			buf=0.0;
		} else {
			buf=2.0*charge[ind0]*fk;
		}	
	
		Vector dd(buf*distABinvmod[k][ind1]*distAB[k][ind1]);
		//derivatives on couple of atoms belonging to groupA and B
		omp_deriv[ind1]+=dd;
		omp_deriv[k]-=dd;
		}
    }

    #pragma omp critical
    for(unsigned i=0; i<getPositions().size(); i++) deriv[i]+=omp_deriv[i];
    #pragma omp barrier

  }

  for(unsigned i=0; i<deriv.size(); ++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (totcharge);
  setBoxDerivatives  (virial);

}
}
}
