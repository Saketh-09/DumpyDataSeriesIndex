class ConversionUtil:
    @staticmethod
    def str_to_int(value, default_value=0):
        try:
            return int(value)
        except ValueError:
            return default_value

    @staticmethod
    def str_to_float(value, default_value=0.0):
        try:
            return float(value)
        except ValueError:
            return default_value

    @staticmethod
    def int_to_str(value):
        return str(value)

    @staticmethod
    def float_to_str(value):
        return str(value)

    @staticmethod
    def bool_to_str(value):
        return "true" if value else "false"

    @staticmethod
    def str_to_bool(value):
        return value.lower() in ["true", "yes", "on", "1"]

    @staticmethod
    def bool_to_int(value):
        return 1 if value else 0

    @staticmethod
    def int_to_bool(value):
        return value != 0

    @staticmethod
    def float_to_int(value):
        return int(value)

    @staticmethod
    def int_to_float(value):
        return float(value)

    @staticmethod
    def float_to_bool(value):
        return bool(value)

    @staticmethod
    def bool_to_float(value):
        return float(value)
