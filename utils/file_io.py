from itertools import islice


def read_data_file_first_line(file_name: str) -> str:
    with open(file_name, "r", encoding="utf-8-sig") as f:
        return f.readline()


def read_data_file_first_five_lines(file_name: str) -> str:
    with open(file_name, "r", encoding="utf-8-sig") as f:
        return " ".join(line.rstrip("\r\n") for line in islice(f, 5))


def read_data_file_all_lines(file_name: str) -> list[str]:
    with open(file_name, "r", encoding="utf-8-sig") as f:
        return f.read().splitlines()
