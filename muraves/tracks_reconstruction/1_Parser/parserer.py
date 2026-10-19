import logging
import tempfile
from typing import NamedTuple, List
from collections import Counter



def parser_new(fileoutp, eventsTot, run, subrun, numEVT, numSk, color):
    eventsTot = int(eventsTot)
    run = int(run) - 1
    subrun = int(subrun) - 1
    
    
    
    #eventsTot = int(infoCounter[1].replace('evts', ''))

    #run = int(infoCounter[2].replace('run', '')) - 1

    #subrun = int(infoCounter[3].replace('sr', '')) - 1

    

    triggerTot = 0 + subrun*int(numEVT) + run*eventsTot
    #brokenF = open('BROKENFILELIST.txt', 'a')

    for ii in range (0,1):
        #fileoutp = sys.argv[1]
        inputF = open(fileoutp, 'r')
        s = inputF.readline()
        linecounter = 0
        CK = 0
        rr = []
        while s:
            s = inputF.readline()
            linecounter = linecounter + 1
        inputF.close()
        rr = s.split('\t')
        if (len(rr)>1) or (len(rr)==0) :
            CK = 1
        else :
            CK = 0
        logging.debug("\n-----" + str(linecounter) + "-----\n")
        lines = open(fileoutp, 'r').readlines()
        rr = lines[160000].split('\t')
        if (rr[0] == '15') :
            if ((linecounter-2) > 160000):
                for i in range(160001, linecounter-1):
                    del lines[160001]
                open(fileoutp, 'w').writelines(lines)
        else :
            CK = 1
        #The following line was commented in the parsing of the NERO, why?? This is the only difference between the three colors.
        #if (CK == 0) and  ((linecounter-2)%int(numSk) == 0) :
        #import pdb; pdb.set_trace()
        if color=='NERO':
            condition = (CK == 0)
        elif color == 'ROSSO' or color == 'BLU' :
            condition = (CK == 0) and  ((linecounter-2)%int(numSk) == 0)

        if condition == True:
            linecounter =160002
            inputF = open(fileoutp, 'r')
            ss = inputF.readlines()
            inputF.close()
            i = 1
            while (i < linecounter-1):
                tss_chk = 0
                rr = ss[i].split('\t')
                for j in range(1,int(numSk)) :
                    if (rr[40] in ss[i+j]) :
                        tss_chk = tss_chk + 1
                if (tss_chk == int(numSk)-1) :
                    i = i + int(numSk)
                else :
                    for j in range(i+tss_chk+1, linecounter) :
                        ss[j-tss_chk-1] = ss[j]
                    linecounter = linecounter - tss_chk - 1
                    if not('2a' in rr[0]) :
                        i = i + int(numSk)
                    else :
                        temp_even = ['','','','','','','','','','','']
                        temp_odd = ['','','','','','','','','','','']
                        for j in range(2, 2*int(numSk), 2):
                            temp_even[int((j-2)/2)] = ss[i+j-1]
                            temp_odd[int((j-2)/2)] = ss[i+j]
                        for j in range(0, int(numSk)-1):
                            ss[i+j+1] = temp_odd[j]
                        for j in range(0,int(numSk)-1):
                            ss[i+j+int(numSk)] = temp_even[j]
                        i = i + int(numSk)*2
            #import pdb; pdb.set_trace()
            with tempfile.NamedTemporaryFile(mode='w+', dir="/tmp", delete=True) as temp_file:
                temp_file.writelines(ss[:linecounter])
                temp_file.seek(0)
                temp_file.readline()  # skip header line
                #import pdb; pdb.set_trace()
            #inputF = open(fileoutp, 'w')
            #for i in range(0, linecounter):
            #    inputF.write(ss[i])
            #inputF.close()

            #mod = fileoutp.replace('slave','mod')
            #mod = sys.argv[2]
            #inputF = open(fileoutp, 'r')
                
                #print(f'Processing file: {mod}')
                offset = 0
                trigger = 0
                # s = inputF.readline()
                p=['']* (int(numSk) * 40)
                #for i in range(int(numSk)*40):
                #    p.append('')
                expected_triggers = (linecounter - 2) // int(numSk)
                #print(f'Reading timestamp to skip', temp_file.readline())
                with tempfile.NamedTemporaryFile(mode='w+', dir="/tmp", delete=False) as temp_outputfile:
                    for s in temp_file:
                        r = s.split('\t')
                        #import pdb; pdb.set_trace()
                        # --- parsing / conversion logic ---
                        h = int('0x'+r[1][2:], 16)

                        p[0+39*offset] = str(h)
                        p[1+39*offset] = r[0]                   
                        h = int('0x'+r[2], 16)
                        p[2+39*offset] = str(h)
                        for i in range (3, 34):
                            if (i <= 17):
                                h = int('0x'+r[i][1:], 16)
                                p[i+39*offset] = str(h)
                            else :
                                h = int('0x'+r[i][2:], 16)
                                p[i+39*offset] = str(h)

                        p[34+39*offset] = r[40].replace('\n', '')
                        h = int('0x'+r[36], 16)
                        p[35+39*offset] = str(h)
                        h = int('0x'+r[37], 16)
                        p[36+39*offset] = str(h)
                        h = int('0x'+r[38][1:], 16)
                        p[37+39*offset] = str(h)
                        h1 = int('0x'+r[34], 16)
                        h2 = int('0x'+r[35], 16)
                        w1 = '{0:b}'.format(h1)
                        w2 = '{0:b}'.format(h2)
                        wh1 = ''
                        wh2 = ''

                        for i in range(16-len(w1)):
                            wh1 = wh1 + '0'
                        for i in range(16-len(w2)):
                            wh2 = wh2 + '0'
                        w1 = wh1 + w1
                        w2 = wh2 + w2
                        tm = '_'
                        for i in range(len(w1)):
                            if(w1[i] == '1'):
                                tm = tm + str(31 - i) + '_'
                        for i in range(len(w2)):
                            if(w2[i] == '1'):
                                tm = tm + str(15 - i) + '_'

                        p[38+39*offset] = tm

                        if (offset < int(numSk)-1) :
                            offset += 1
                        else :
                            temp_outputfile.write(str(triggerTot) + '\t')
                            temp_outputfile.write('\t'.join(p[:int(numSk)*39]) + '\n')
                            offset = 0
                            trigger += 1
                            triggerTot += 1

                        
                        if trigger >= expected_triggers:
                            break
                    #if (trigger==(linecounter-2)/int(numSk)) :
                    #    s = inputF.readline()
                    #    s = inputF.readline()
                #import pdb; pdb.set_trace()
                #inputF.close()
                #outputF.close()

        else :
            print("\n STATUS : Skip File\n")
            
            #brokenF.write(fileoutp + '\n')
    return CK, temp_outputfile.name
    #brokenF.close()
def single_bit_flip(a: int, b: int) -> bool:
    """
    Returns True if a and b differ by exactly one bit.
    """
    return (a ^ b).bit_count() == 1




class ValidationResult(NamedTuple):
    """Result of event number spacing validation."""
    errors: List[tuple[int, str]]  # (line_number, error_message)
    has_mismatches: bool
    has_unrecoverable_errors: bool
    mismatches_counter: int
    unrecoverable_mismatches_counter: int
    mismatches_per_line: dict # line_number -> count of mismatches on that line
    unrecoverable_mismatches_details: dict 

def insertion_test(a: int, b: int) -> bool:
    """
    Returns True if inserting exactly one bit (0 or 1)
    anywhere into the binary string of 'a' can make it equal to 'b',
    or vice versa.
    """

    ba = bin(a)[2:]
    bb = bin(b)[2:]

    # One must be exactly one bit longer than the other
    if abs(len(ba) - len(bb)) != 1:
        return False

    # Identify shorter and longer
    short, long = (ba, bb) if len(ba) < len(bb) else (bb, ba)

    # Try inserting '0' or '1' into every possible position of short
    for i in range(len(short) + 1):
        for bit in ('0', '1'):
            new = short[:i] + bit + short[i:]
            if new == long:
                return True

    return False

def check_subevent_number_spacing(file_obj, block_size=39, line_to_print=100) -> ValidationResult:
    """
    Validates that event numbers repeat every `block_size` tokens in a file object.
    
    For each row:
      - Extracts the event number from the second field.
      - Verifies that the same event number appears at token positions: 1, 1+block_size, 1+2*block_size, ...

    Args:
        file_obj: File-like object to validate.
        block_size: Expected spacing between event numbers (default: 39).
        line_to_print: Line number for debug logging (default: 100).

    Returns:
        ValidationResult containing:
            - errors: List of (line_number, error_message) tuples.
            - has_mismatches: True if any mismatch found.
            - has_unrecoverable_errors: True if mismatch failed the insertion_test (real error).
    """
    logging.debug("Running event number spacing check...")
    errors = []
    has_mismatches = False
    has_unrecoverable = False
    mismatches_counter = 0
    unrecoverable_mismatches_counter = 0
    mismatches_per_line = {
        "line_number": [],
        "event_number": [],
        "subevent_number_expected":[],
        "subevent_number_found":[],
        "block_idx": [],
        "bit_flip": [],
        "bit_missing": [],
    }
    unrecoverable_mismatches_details = {
        "line_number": [],
        "event_number": [],
        "block_index": [],
        "subevent_number_expected": [],
        "subevent_number_found": []
    }  # List of (line_number, block_index, expected, found)

    #previous_subevent_number = 0
    for line_number, line in enumerate(file_obj, start=1):
        #mismatches_counter_per_line = 0
        tokens = line.strip().split("\t")
        if len(tokens) < 2:
            errors.append((line_number, "Row too short to contain event number"))
            continue
        try:
            event_number = int(tokens[0])
            # take the number that is the most frequent among the the i*block_size with i from 1 to block_size
            subevent_number = Counter([int(tokens[i]) for i in range(1, len(tokens), block_size)]).most_common(1)[0][0]
            #logging.debug(f'Event number: {subevent_number}')
            #subevent_number = int(tokens[1])
            if line_number == line_to_print:
                logging.debug(f'Event number: {subevent_number}')
        except ValueError:
            errors.append((line_number, "Invalid event number in column 2"))
            continue
        
        
        # Track if sequential number was corrected
        # Vorrei cambiare questo scegliendo come subevent_number quello più frequente tra i blocks
        #if line_number > 1:
        #    current_subevent_number_expected = previous_subevent_number + 1
        #    if current_subevent_number_expected == 4096: # wrap around case
        #        current_subevent_number_expected = 0
        #    if line_number == line_to_print:
        #        logging.debug(f'Previous subevent number: {previous_subevent_number}')
        #        logging.debug(f'Expected subevent number: {current_subevent_number_expected}')
        #    if int(subevent_number) != int(current_subevent_number_expected):
        #        logging.debug(f"Mismatch found at line {line_number}: expected {current_subevent_number_expected}, got {subevent_number}, setting right eventnumber!")
        #        subevent_number = current_subevent_number_expected  # Sostituisci con il valore atteso
        #        mismatches_counter += 1
        #        mismatches_counter_per_line += 1
        #        has_mismatches = True
        
        # Check event number repeats at block_size intervals: always start from second block
        # since subevent_number is always taken from first field (position 1)
        idx = 1 
        block_index = 1
        while idx < len(tokens):
            try:
                value = int(tokens[idx])
                if line_number == line_to_print:
                    logging.debug(value)
                    logging.debug(subevent_number)
            except ValueError:
                errors.append((line_number, f"Non-integer token at position {idx}"))
                break
            
            if value != subevent_number:   
                has_mismatches = True
                mismatches_counter += 1
                #mismatches_counter_per_line += 1
                mismatches_per_line["line_number"].append(line_number)
                mismatches_per_line["event_number"].append(event_number)
                mismatches_per_line["subevent_number_expected"].append(subevent_number)
                mismatches_per_line["subevent_number_found"].append(value)
                mismatches_per_line["block_idx"].append(block_index)
                #mismatches_per_line["number_mismatches"].append(mismatches_counter_per_line)

                #block_mismatch = block_index
                # Check if this could be a single-bit insertion error
                if not single_bit_flip(value, subevent_number):
                    if not insertion_test(value, subevent_number):
                        has_unrecoverable = True
                        unrecoverable_mismatches_counter += 1
                        unrecoverable_mismatches_details["line_number"].append(line_number)
                        unrecoverable_mismatches_details["event_number"].append(event_number)
                        unrecoverable_mismatches_details["block_index"].append(block_index)
                        unrecoverable_mismatches_details["subevent_number_expected"].append(subevent_number)
                        unrecoverable_mismatches_details["subevent_number_found"].append(value)
                        #logging.debug(f"Unrecoverable error for line {line_number}")
                        errors.append(
                            (f"Unrecoverable mismatch at block {block_index} of line {line_number}: "
                             f"expected subevent number {subevent_number}, but found {value}")
                        )
                    else:
                        mismatches_per_line["bit_missing"].append(True)
                        mismatches_per_line["bit_flip"].append(False)
                else:
                    mismatches_per_line["bit_missing"].append(False)
                    mismatches_per_line["bit_flip"].append(True)
                # else: mismatch passes insertion_test, so it's "recoverable" (don't record)
                #break
            
            idx += block_size
            block_index += 1

        #if mismatches_counter_per_line > 0:


        #previous_subevent_number = subevent_number
        
        

            
    logging.debug(f"Finished check. Found {len(errors)} errors.")
    return ValidationResult(errors=errors, has_mismatches=has_mismatches, has_unrecoverable_errors=has_unrecoverable, mismatches_counter=mismatches_counter, unrecoverable_mismatches_counter=unrecoverable_mismatches_counter, mismatches_per_line=mismatches_per_line, unrecoverable_mismatches_details=unrecoverable_mismatches_details)



