import pprint

grades = {
    "Biology" : 4,
    "Chemistry" : 2,
    "Maths" : 3,
    "Physics" : 3,
    "Geography" : 5,
}

pprint.pp((pprint.pformat(grades,compact=False)))