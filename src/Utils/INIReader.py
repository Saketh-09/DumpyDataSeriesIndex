class INIReader:
    def __init__(self, filename=None, buffer=None):
        self._values = {}
        self._error = 0
        if filename:
            self._error = self._parse_file(filename)
        elif buffer:
            self._error = self._parse_string(buffer)

    def _parse_file(self, filename):
        try:
            with open(filename, 'r') as file:
                content = file.read()
                return self._parse_string(content)
        except FileNotFoundError:
            return 1  # Error code for file not found

    def _parse_string(self, content):
        lines = content.split('\n')
        current_section = None
        for line in lines:
            line = line.strip()
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1]
            elif '=' in line:
                key, value = line.split('=', 1)
                key = key.strip()
                value = value.strip()
                if current_section:
                    key = f'{current_section}={key}'
                self._values[key] = value
        return 0  # No error

    def parse_error(self):
        return self._error

    def get(self, section, name, default_value=None):
        key = f'{section}={name}'
        return self._values.get(key, default_value)

    def get_string(self, section, name, default_value=None):
        value = self.get(section, name, default_value)
        return str(value) if value is not None else default_value

    def get_integer(self, section, name, default_value=0):
        value = self.get(section, name)
        return int(value) if value is not None else default_value

    def get_real(self, section, name, default_value=0.0):
        value = self.get(section, name)
        return float(value) if value is not None else default_value

    def get_boolean(self, section, name, default_value=False):
        value = self.get(section, name)
        if value is not None:
            value = value.lower()
            if value in ['true', 'yes', 'on', '1']:
                return True
            elif value in ['false', 'no', 'off', '0']:
                return False
        return default_value

    def has_section(self, section):
        return any(key.startswith(f'{section}=') for key in self._values)

    def has_value(self, section, name):
        key = f'{section}={name}'
        return key in self._values

    def _value_handler(self, section, name, value):
        key = f'{section}={name}'
        if key in self._values and self._values[key]:
            self._values[key] += '\n'
        self._values[key] += value if value else ''

    def load_from_buffer(self, buffer):
        return self._parse_string(buffer)


