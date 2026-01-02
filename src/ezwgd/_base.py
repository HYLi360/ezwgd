""""""

import platform

version = "0.1.0"
nickname = "Aegis"


def base_env():
    print(f"""System          {platform.platform()}
Python Path     {platform.python_version()}
Python Version  {platform.python_version()}""")
