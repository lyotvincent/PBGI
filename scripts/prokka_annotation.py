import subprocess


def annotation(input_file, result_dir, conf):
    print("Begin Prokka")
    print("input_file="+input_file)
    prokka_conf = conf['denovo']['prokka']
    command_line = 'prokka --outdir '+result_dir+'/denovo/prokka_annotation '
    if not prokka_conf['--addgenes'] == 0:
        command_line += '--addgenes '
    if not prokka_conf['--addmrna'] == 0:
        command_line += '--addmrna '
    if not prokka_conf['--evalue'] == 0:
        command_line += '--evalue '+str(prokka_conf['--evalue'])+' '
    if not prokka_conf['--coverage'] == 0:
        command_line += '--coverage '+str(prokka_conf['--coverage'])+' '
    if not prokka_conf['--cpus'] == 0:
        command_line += '--cpus '+str(prokka_conf['--cpus'])+' '
    if not prokka_conf['--mincontiglen'] == 0:
        command_line += '--mincontiglen '+str(prokka_conf['--mincontiglen'])+' '
    if not prokka_conf['--rfam'] == 0:
        command_line += '--rfam '
    if not prokka_conf['--norrna'] == 0:
        command_line += '--norrna '
    if not prokka_conf['--notrna'] == 0:
        command_line += '--notrna '
    if not prokka_conf['--rnammer'] == 0:
        command_line += '--rnammer '
    command_line += input_file
    subprocess.run(command_line, shell=True, check=True)
    print("End Prokka")
