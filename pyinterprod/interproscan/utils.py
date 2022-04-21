def int_to_upi(i):
    return f"UPI{i:010x}".upper()


def upi_to_int(upi):
    return int(upi[3:], 16)


def range_jobs(from_upi: str, to_upi: str, step: int):
    start = upi_to_int(from_upi)
    stop = upi_to_int(to_upi) + 1
    for i in range(start, stop, step):
        yield int_to_upi(i), int_to_upi(i + step - 1)
