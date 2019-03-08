import sys
import re


# Iterates through stdin looking for the string that matches kraken2's
# classification report. Output will look like this:
#   Loading database information... done.
#   1452549 sequences (680.14 Mbp) processed in 13.338s (6534.0 Kseq/m, 3059.49 Mbp/m).
#       1444500 sequences classified (99.45%) <---- matches this line
#       8049 sequences unclassified (0.55%)
#
# Then finds the percentage classified, and returns 0 only if the percentage
# classfied is greater than or equal to 90%.
def main():
    classified = re.compile(r'(?<!(un))classified \([0-9]+\.[0-9]*%\)')
    for line in sys.stdin:
        match = classified.search(line)
        if match:
            perc = re.compile(r'[0-9]+\.[0-9]*')
            percent = perc.search(match.group(0))
            tbperc = float(percent.group(0)) / 100

            if tbperc < 0.9:
                exit(1)

            break


if __name__ == '__main__':
    main()
