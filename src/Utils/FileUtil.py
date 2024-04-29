import os

class FileUtil:
    @staticmethod
    def file_exists(filename):
        return os.path.isfile(filename)

    @staticmethod
    def directory_exists(directory):
        return os.path.isdir(directory)

    @staticmethod
    def create_directory(directory):
        try:
            os.makedirs(directory)
            return True
        except OSError:
            return False

    @staticmethod
    def remove_file(filename):
        try:
            os.remove(filename)
            return True
        except OSError:
            return False

    @staticmethod
    def rename_file(old_filename, new_filename):
        try:
            os.rename(old_filename, new_filename)
            return True
        except OSError:
            return False

    @staticmethod
    def read_file(filename):
        try:
            with open(filename, 'r') as file:
                return file.read()
        except FileNotFoundError:
            return None

    @staticmethod
    def write_file(filename, content):
        try:
            with open(filename, 'w') as file:
                file.write(content)
            return True
        except OSError:
            return False
