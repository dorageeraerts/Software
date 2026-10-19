
__all__= ['parse_runs', 'batch_runs', 'parse_batch_indices']

# ******************** Funzioni helper  ********************
#Function to convert runs of config file to list of runs
def parse_runs(run_str):
    runs = []
    try:
        for part in run_str.split(','):
            part = part.strip()
            if '-' in part:
                start, end = map(int, part.split('-'))
                runs.extend(range(start, end + 1))
            else:
                runs.append(int(part))
        
    except:
        runs.append(run_str)
    
    return [str(r) for r in runs] 

def batch_runs(runs, n=5):
    """Split list of runs into chunks of size n"""
    for i in range(0, len(runs), n):
        yield runs[i:i+n]

def parse_batch_indices(batch_string):
    batch_list = []
    parts = batch_string.split(',')

    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            batch_list.extend(range(start, end + 1))
        else:
            batch_list.append(int(part))

    return sorted(set(batch_list))  # Rimuove duplicati e ordina