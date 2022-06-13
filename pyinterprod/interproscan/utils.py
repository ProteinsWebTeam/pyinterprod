def int_to_upi(i: int, digits: int = 10):
    return f"UPI{i:0{digits}x}".upper()


def upi_to_int(upi):
    return int(upi[3:], 16)


def range_upi(from_upi: str, to_upi: str, step: int):
    start = upi_to_int(from_upi)
    stop = upi_to_int(to_upi) + 1
    digits = len(from_upi) - 3  # number of digits after "UPI"
    for i in range(start, stop, step):
        yield (
            int_to_upi(i, digits=digits),
            min(to_upi, int_to_upi(i + step - 1, digits=digits))
        )
