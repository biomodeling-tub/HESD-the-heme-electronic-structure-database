import os
import sys

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} SEARCH_TERM")
    sys.exit(1)
search_term = sys.argv[1]
count = sum(
    search_term in f
    for f in os.listdir('.')
    if os.path.isfile(f)
)
print(f"Number of files containing '{search_term}' in the name: {count}")