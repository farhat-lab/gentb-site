import sys


# check what proportion of depths are >= 10
# want at least 95% >= 10
def main():
    num_lines = 0
    valid = 0

    for line in sys.stdin:
        num_lines += 1

        # samtools depth -a format has depth in third column
        depth = int(line.split('\t')[2])
        if depth >= 10:
            valid += 1

    # nonzero exit for GenTB Pipeline
    if (valid / num_lines) < 0.95:
        exit(1)


if __name__ == '__main__':
    main()
