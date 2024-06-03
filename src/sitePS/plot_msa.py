from pymsaviz import MsaViz
import os


def msa_plot(out_csv,
             png,
             namelist: str = 'Basic_101species_name.list',
             plot_length: int = 120,
             plot_dpi: int = 500):
    try:
        with open(namelist) as f:
            file = f.readlines()
        dict_name = {}
        file.pop(0)
        for line in file:
            line = line.rstrip('\n')
            uid, v, k = line.split(',')
            dict_name[k] = v
    except Exception as e:
        print(e)
        print('\nMaybe there is no Basic_101species_name.list')
        return
    with open(out_csv) as f:
        file = f.readlines()
    file.pop(0)
    sites = file.pop().split(',')[2:]
    length = len(sites)
    with open(out_csv + '.fas', 'w') as w:
        for line in file:
            line = line.rstrip('\n')
            name = line.split(',')[1]
            seq = ''.join(line.split(',')[2:])
            try:
                PS = dict_name[name]
                name = name.replace(' ', '_')
                w.write('>' + name + '-' + PS + '\n')
            except Exception:
                name = name.replace(' ', '_')
                w.write('>' + name + '\n')
            w.write(seq + '\n')

    mv = MsaViz(out_csv + '.fas', wrap_length=plot_length, show_consensus=True)
    for item in range(1, length + 1):
        mv.add_text_annotation((item, item), text=sites[item - 1], text_size=7)
    mv.savefig(png, dpi=plot_dpi)
    os.remove(out_csv + '.fas')



