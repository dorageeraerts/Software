import logging
import tempfile

def parser_nero(fileoutp, mod, numEVT, numSk):


    infoCounter = fileoutp.split('_')
    eventsTot = int(infoCounter[1].replace('evts', ''))
    run = int(infoCounter[2].replace('run', '')) - 1
    subrun = int(infoCounter[3].replace('sr', '')) - 1
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
        print("\n-----" + str(linecounter) + "-----\n")
        lines = open(fileoutp, 'r').readlines()
        rr = lines[160000].split('\t')
        if (rr[0] == '15') :
            if ((linecounter-2) > 160000):
                for i in range(160001, linecounter-1):
                    del lines[160001]
                open(fileoutp, 'w').writelines(lines)
        else :
            CK = 1
        #if (CK == 0) and  ((linecounter-2)%int(numSk) == 0) :
        if (CK == 0) :
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
            
            inputF = open(fileoutp, 'w')
            for i in range(0, linecounter):
                inputF.write(ss[i])
            inputF.close()

            #mod = fileoutp.replace('slave','mod')
            #mod = sys.argv[2]
            inputF = open(fileoutp, 'r')
            outputF = open(mod, 'a')
            offset = 0
            trigger = 0
            s = inputF.readline()
            p=['']
            for i in range(int(numSk)*40):
                p.append('')

            while s:
                s = inputF.readline()
                r = s.split('\t')
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
                    offset = offset + 1
                else :
                    outputF.write(str(triggerTot) + '\t')
                    for i in range(int(numSk)*39-1):
                        outputF.write(p[i] + '\t')
                    outputF.write(p[int(numSk)*39-1] + '\n')
                    offset = 0
                    trigger = trigger + 1
                    triggerTot = triggerTot + 1

                if (trigger==(linecounter-2)/int(numSk)) :
                    s = inputF.readline()
                    s = inputF.readline()

            inputF.close()
            outputF.close()


        else :
            print("\n STATUS : Skip File\n")
            
            #brokenF.write(fileoutp + '\n')

    return CK


def parser_nero_new(fileoutp, mod, numEVT, numSk):


    infoCounter = fileoutp.split('_')

    eventsTot = int(infoCounter[1].replace('evts', ''))

    run = int(infoCounter[2].replace('run', '')) - 1

    subrun = int(infoCounter[3].replace('sr', '')) - 1

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
        #if (CK == 0) and  ((linecounter-2)%int(numSk) == 0) :
        if (CK == 0) :
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
            with tempfile.NamedTemporaryFile(mode='w+', delete=True) as temp_file:
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
                with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_outputfile:
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


def parser_rosso(fileoutp, mod, numEVT, numSk):
    infoCounter = fileoutp.split('_')
    eventsTot = int(infoCounter[1].replace('evts', ''))
    run = int(infoCounter[2].replace('run', '')) - 1
    subrun = int(infoCounter[3].replace('sr', '')) - 1
    triggerTot = 0 + subrun*int(numEVT) + run*eventsTot
    #brokenF = open('BROKENFILELIST.txt', 'a')
    
    for ii in range (0,1):
        inputF = open(fileoutp, 'r')
        s = inputF.readline()
        linecounter = 0
        CK = 0
        rr = []
        while s:
            s = inputF.readline()
            linecounter = linecounter + 1
        rr = s.split('\t')
        if (len(rr)>1) or (len(rr)==0) :
            CK = 1
        else :
            CK = 0
        inputF.close()
        print("\n-----" + str(linecounter) + "-----\n")
        if (CK == 0) and  ((linecounter-2)%int(numSk) == 0) :
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
    
            inputF = open(fileoutp, 'w')
            for i in range(0, linecounter):
                inputF.write(ss[i])
            inputF.close()
    
            #mod = fileoutp.replace('slave','mod')
            inputF = open(fileoutp, 'r')
            outputF = open(mod, 'a')
            offset = 0
            trigger = 0
            s = inputF.readline()
            p=['']
            for i in range(int(numSk)*40):
                p.append('')
    
            while s:
                s = inputF.readline()
                r = s.split('\t')
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
                    offset = offset + 1
                else :
                    outputF.write(str(triggerTot) + '\t')
                    for i in range(int(numSk)*39-1):
                        outputF.write(p[i] + '\t')
                    outputF.write(p[int(numSk)*39-1] + '\n')
                    offset = 0
                    trigger = trigger + 1
                    triggerTot = triggerTot + 1
    
                if (trigger==(linecounter-2)/int(numSk)) :
                    s = inputF.readline()
                    s = inputF.readline()
    
            inputF.close()
            outputF.close()
    
        else :
            print("\n STATUS : Skip File\n")
            #brokenF.write(fileoutp + '\n')
    
    #brokenF.close()
    return CK

def parser_blu(fileoutp, mod, numEVT, numSk):


    infoCounter = fileoutp.split('_')
    eventsTot = int(infoCounter[1].replace('evts', ''))
    run = int(infoCounter[2].replace('run', '')) - 1
    subrun = int(infoCounter[3].replace('sr', '')) - 1
    triggerTot = 0 + subrun*int(numEVT) + run*eventsTot
    #brokenF = open('BROKENFILELIST.txt', 'a')
    
    for ii in range (0,1):
        fileoutp = fileoutp
        inputF = open(fileoutp, 'r')
        s = inputF.readline()
        linecounter = 0
        CK = 0
        rr = []
        while s:
            s = inputF.readline()
            linecounter = linecounter + 1
        rr = s.split('\t')
        if (len(rr)>1) or (len(rr)==0) :
            CK = 1
        else :
            CK = 0
        inputF.close()
        print("\n-----" + str(linecounter) + "-----\n")
        if (CK == 0) and  ((linecounter-2)%int(numSk) == 0) :
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
    
            inputF = open(fileoutp, 'w')
            for i in range(0, linecounter):
                inputF.write(ss[i])
            inputF.close()
    
            #mod = fileoutp.replace('slave','mod')
            inputF = open(fileoutp, 'r')
            outputF = open(mod, 'a')
            offset = 0
            trigger = 0
            s = inputF.readline()
            p=['']
            for i in range(int(numSk)*40):
                p.append('')
    
            while s:
                s = inputF.readline()
                r = s.split('\t')
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
                    offset = offset + 1
                else :
                    outputF.write(str(triggerTot) + '\t')
                    for i in range(int(numSk)*39-1):
                        outputF.write(p[i] + '\t')
                    outputF.write(p[int(numSk)*39-1] + '\n')
                    offset = 0
                    trigger = trigger + 1
                    triggerTot = triggerTot + 1
    
                if (trigger==(linecounter-2)/int(numSk)) :
                    s = inputF.readline()
                    s = inputF.readline()
    
            inputF.close()
            outputF.close()
    
        else :
            print("\n STATUS : Skip File\n")
            #brokenF.write(fileoutp + '\n')
    
    #brokenF.close()
    return CK
    




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

def check_event_number_spacing(filename, block_size=39, line_to_print =1):
    """
    For each row:
      - Extract the event number from the second field.
      - Verify that the same event number appears again every `block_size` items.

    A row is considered valid if, for k = 1, 2, ..., it holds that:
        tokens[k * block_size + 1] == event_number

    Returns a list of (line_number, message) errors.
    """
    logging.debug("Runing relative events number check...")
    errors = []

    with open(filename, "r") as f:
        for line_number, line in enumerate(f, start=1):

            tokens = line.strip().split("\t")

            if len(tokens) < 2:
                errors.append((line_number, "Row too short to contain event number"))
                continue

            try:

                event_number = int(tokens[1])
                if line_number==line_to_print:
                    print('Event number: ', event_number)
            except ValueError:
                errors.append((line_number, "Invalid event number in column 2"))
                continue

            # Now check repeating event numbers every block_size columns
            # The event number should appear at token positions: 1, 1+block_size, 1+2*block_size, ...
            idx = 1
            block_index = 1

            while idx < len(tokens):
                try:
                    value = int(tokens[idx])
                    if line_number==line_to_print:
                        print(value)
                except ValueError:
                    errors.append((line_number, f"Non-integer token at position {idx}"))
                    break

                if value != event_number and not insertion_test:
                    errors.append(
                        (line_number,
                         f"Event number mismatch at block {block_index}: "
                         f"found {value}, expected {event_number}")
                    )
                    break

                idx += block_size
                block_index += 1

    return errors

