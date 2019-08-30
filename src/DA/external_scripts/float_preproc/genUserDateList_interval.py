from datetime import date, datetime, timedelta

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

for result in perdelta(date(2015, 01, 3), date(2015, 12, 31), timedelta(days=3)):
    print result
