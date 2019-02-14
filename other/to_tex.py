def main():
    cols = 7
    i = 0
    firstLine = True
    res = '\\begin{tabular}{|' + ('c|' * cols) + '}\n'
    res += '\hline\n'
    for l in open('in.txt', 'r').readlines():
        l = l.strip().replace('±', '\pm').replace('–', '-')
        l = '$' + l + '$'
        if i < cols - 1:
            res += l + ' & '
            i += 1
        else:
            res += l + ' \\\\\n'
            i = 0
            if firstLine:
                res += '\hline\n'
                firstLine = False

    res += '\\hline\n'
    res += '\end{tabular}'
    open('table4.txt', 'w').write(res)
    print(res)


if __name__ == '__main__':
    main()