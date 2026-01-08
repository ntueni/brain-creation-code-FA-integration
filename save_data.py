# Globale Variable, um das Array zu speichern
_saved_data = []

def save_data(data):
    """Speichert das Array in einer globalen Variablen."""
    global _saved_data
    _saved_data = data

def get_data():
    """Gibt das gespeicherte Array zurück."""
    return _saved_data