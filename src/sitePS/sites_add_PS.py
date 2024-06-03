def site_add_PS(sites_csv, output):
    with open(sites_csv) as f:
        file = f.readlines()
    bak = file.copy()
    file.pop(0)
    human = file.pop(0)
    human_seq = human.rstrip('\n').split(',')[2:]
    human_stratum = human.split(',')[0]
    stratum = {}
    stratum_list = []
    total_other_seq = []
    for line in file:
        line = line.rstrip('\n')
        line = line.split(',')
        stratum_list.append(line[0])
        if line[0] not in stratum:
            stratum[line[0]] = [line[1]]
        else:
            stratum[line[0]].append(line[1])
        other_seq = line[2:]
        total_other_seq.append(other_seq)
    stratum_list.append(human_stratum)
    no_repeat_stream_list = sorted(list(set(stratum_list)), key=lambda x: int(x), reverse=True)
    record = []
    for i in range(len(human_seq)):
        aa_human = human_seq[i]
        other_aa = [item[i] for item in total_other_seq]
        length_record = len(record)
        for aa_idx in range(len(other_aa)):
            if other_aa[aa_idx] == aa_human:
                continue
            else:
                current_stratum = stratum_list[aa_idx]
            if len(stratum[current_stratum]) == 1:
                record.append(no_repeat_stream_list[no_repeat_stream_list.index(current_stratum) - 1])
                break
            else:
                start = stratum_list.index(current_stratum)
                if aa_human in other_aa[start: start + len(stratum[current_stratum])]:
                    continue
                else:
                    record.append(no_repeat_stream_list[no_repeat_stream_list.index(current_stratum) - 1])
                    break
        if len(record) == length_record:
            if other_aa[0] != aa_human:
                record.append(no_repeat_stream_list[0])  # human PS
            else:
                record.append(no_repeat_stream_list[-1])


    xx = [',' + item for item in record]
    with open(output, 'w') as w:
        w.writelines(bak)
        w.write(',')
        w.writelines(xx)
        w.write('\n')

