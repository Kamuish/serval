

def create_print_file(filename, *args):
    with open(filename + options[chosen_option], mode = 'w') as file:
        print(*args, file = file)


if __name__ == '__main__':
    create_print_file('teste.txt', [''],0, 'a','b')