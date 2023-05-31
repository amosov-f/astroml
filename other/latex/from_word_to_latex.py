import re

def main():
    with open('8.txt') as f:
        l = 0
        for line in f.readlines():
            line = line.strip()
            if line:
                line = line.replace(' ± ', '±')
                parts = re.split('\s+', line)
                if l == 0:
                    print('|c' * len(parts) + '|')
                    print('\hline')
                for i in range(len(parts)):
                    part = parts[i]
                    if part.startswith('s') or part.startswith('t'):
                        part = f'{part[0]}_{{{part[1:]}}}'
                    if part.startswith('b'):
                        part = f'\\textbf{{{part[1:]}}}'
                    print(f'${part}$', end='')
                    if i < len(parts) - 1:
                        print(' & ', end='')
                print(' \\\\')
                print('\hline')
                l += 1

if __name__ == '__main__':
    main()