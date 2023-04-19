
import socket
import os

def autodetect():

    """
    Returns
    -------
    bool
    	True if current platform matches, otherwise False
    """

    if socket.gethostname() == 'f-dahu':
        return True

    return False


if __name__ == "__main__":
    print("Autodetect: "+str(autodetect()))

