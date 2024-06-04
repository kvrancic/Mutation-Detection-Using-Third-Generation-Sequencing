import csv

def read_results(filename):
    results = []
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            results.append((row[0], int(row[1]), row[2]))
    return results

def evaluate_results(my_results, ref_results, log_filename):
    total_points = 0
    max_points = 0

    x_points = 1
    y_points = 1
    z_points = 1

    with open(log_filename, 'w') as log_file:
        log_file.write("Evaluation Log\n")
        log_file.write("====================\n")

        my_dict = {pos: (type_, base) for type_, pos, base in my_results}

        for ref_type, ref_pos, ref_base in ref_results:
            max_points += x_points + y_points + z_points
            log_file.write(f"Checking reference result: {ref_type},{ref_pos},{ref_base}\n")
            if ref_pos in my_dict:
                my_type, my_base = my_dict[ref_pos]
                log_file.write(f"  Found matching position in my results: {my_type},{ref_pos},{my_base}\n")
                total_points += x_points
                if my_type == ref_type:
                    total_points += y_points
                    if my_base == ref_base:
                        total_points += z_points
                        log_file.write(f"    All match: +{x_points + y_points + z_points} points\n")
                    else:
                        log_file.write(f"    Type match, but base mismatch: +{x_points + y_points} points\n")
                else:
                    log_file.write(f"    Position match, but type mismatch: +{x_points} points\n")
            else:
                log_file.write(f"  No matching position in my results\n")

        log_file.write("====================\n")
        accuracy = (total_points / max_points) * 100
        log_file.write(f"Total points: {total_points}\n")
        log_file.write(f"Max points: {max_points}\n")
        log_file.write(f"Accuracy: {accuracy:.2f}%\n")

    print(f"Accuracy: {accuracy:.2f}%")

# Replace with the actual file paths
my_results_file = "data/mutations.csv"
ref_results_file = "data/lambda_mutated.csv"
log_file = "evaluation_log.txt"

my_results = read_results(my_results_file)
ref_results = read_results(ref_results_file)
evaluate_results(my_results, ref_results, log_file)
