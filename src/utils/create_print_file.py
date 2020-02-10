

def create_print_file(filename, *args):
    with open(filename, mode = 'a') as file:
        print(*args, file = file)


if __name__ == '__main__':
    create_print_file('teste.txt', [''],0, 'a','b')