import os
import sys
if len(sys.argv) < 3:
    print(f"Usage: {sys.argv[0]} SEARCH_TERM EXCLUDE_TERM")
    sys.exit(1)
search_term = sys.argv[1]
exclude_term = sys.argv[2]
count = 0
for f in os.listdir('.'):
    if os.path.isfile(f):
        if search_term in f and exclude_term not in f:
            count += 1
print(f"Files containing '{search_term}' but not '{exclude_term}': {count}")