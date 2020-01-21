class SnakeModule:
    name = ''
    install_dir = ''
    snakefiles = []
    
    def __init__(self, name, install_dir, snakefiles):
        self.name = name
        self.install_dir = install_dir
        self.snakefiles = snakefiles
