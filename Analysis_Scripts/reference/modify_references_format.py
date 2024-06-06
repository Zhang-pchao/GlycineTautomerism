import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
import titlecase
from pyiso4.ltwa import Abbreviate

# Create an abbreviator (using the default LTWA)
abbreviator = Abbreviate.create()

# List of words that should remain in lowercase
lowercase_words = ['from', 'up','into','through', 'than']

def title_case_title(title):
    # Apply titlecase function
    titled = titlecase.titlecase(title)
    # Post-process to ensure certain words remain lowercase
    words = titled.split()
    corrected_title = [words[0]] + [
        word if word.lower() not in lowercase_words else word.lower() for word in words[1:]
    ]
    return ' '.join(corrected_title)

def abbreviate_journal(journal):
    # Check if the journal name is already abbreviated (contains periods)
    if '.' in journal:
        return journal
    abbreviated = abbreviator(journal, remove_part=True)
    return titlecase.titlecase(abbreviated)

def process_bib_file(input_file, output_file, title_log_file, journal_log_file):
    # Read the input .bib file
    with open(input_file, 'r') as bibfile:
        parser = BibTexParser()
        bib_database = bibtexparser.load(bibfile, parser=parser)

    title_changes = []
    journal_changes = []

    # Process each entry
    for entry in bib_database.entries:
        if 'title' in entry:
            original_title = entry['title']
            # Apply title case to the title
            modified_title = title_case_title(original_title)
            if original_title != modified_title:
                title_changes.append(f"Original: {original_title}\nModified: {modified_title}\n")
                entry['title'] = modified_title

        if 'journal' in entry:
            original_journal = entry['journal']
            # Abbreviate the journal name and convert to title case
            modified_journal = abbreviate_journal(original_journal)
            if original_journal != modified_journal:
                journal_changes.append(f"Original: {original_journal}\nModified: {modified_journal}\n")
                entry['journal'] = modified_journal

    # Write the changes to the output .bib file
    writer = BibTexWriter()
    with open(output_file, 'w') as bibfile:
        bibtexparser.dump(bib_database, bibfile, writer=writer)

    # Log the title changes
    with open(title_log_file, 'w') as logfile:
        logfile.write("\n".join(title_changes))

    # Log the journal changes
    with open(journal_log_file, 'w') as logfile:
        logfile.write("\n".join(journal_changes))

# Define file paths
input_bib_file = 'references.bib'
output_bib_file = 'references_changed.bib'
title_log_file = 'change_title.log'
journal_log_file = 'change_journal.log'

# Process the .bib file
process_bib_file(input_bib_file, output_bib_file, title_log_file, journal_log_file)