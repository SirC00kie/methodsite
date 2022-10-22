import json


def write_to_file(str, url):
    file = open(url, 'w')
    file.write(str)
    file.close()


def read_from_file(url):
    file = open(url, 'r')
    text = file.read()
    file.close()
    return text


def get_json(url):
    f = open(url)
    data = json.load(f)
    f.close()
    return data


def set_json(tdict, url):
    with open(url, 'w') as fp:
        json.dump(tdict, fp, indent=4, ensure_ascii=False)


