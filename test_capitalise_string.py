def capitalise_string(s):
    if not isinstance(s, str):
        raise TypeError('Not a string.')
    return s.capitalize()

def test_capitalise_string():
    assert capitalise_string('test')
