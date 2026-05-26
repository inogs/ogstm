from datetime import datetime, timedelta

def create_ave_times(start_str, end_str, step_value, step_unit, output_file="2.aveTimes"):
    fmt = "%Y%m%d-%H:%M:%S"
    start = datetime.strptime(start_str, fmt)
    end = datetime.strptime(end_str, fmt)

    if step_value <= 0:
        raise ValueError("step_value deve essere positivo")

    step_unit = step_unit.lower()

    if step_unit in ["second", "seconds", "sec", "s"]:
        delta = timedelta(seconds=step_value)
    elif step_unit in ["minute", "minutes", "min", "m"]:
        delta = timedelta(minutes=step_value)
    elif step_unit in ["hour", "hours", "h"]:
        delta = timedelta(hours=step_value)
    elif step_unit in ["day", "days", "d"]:
        delta = timedelta(days=step_value)
    else:
        raise ValueError("step_unit deve essere uno tra: seconds, minutes, hours, days")

    current = start

    with open(output_file, "w") as f:
        while current <= end:
            f.write(current.strftime("%Y%m%d-%H:%M:%S") + "\n")
            current += delta

if __name__ == "__main__":
    start_date = "20000101-00:00:00"
    end_date = "20010101-00:00:00"
    step_value = 1
    step_unit = "days"

    create_ave_times(start_date, end_date, step_value, step_unit)
    print("File 2.aveTimes creato con successo.")
