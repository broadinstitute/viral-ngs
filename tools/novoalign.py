# /idi/sabeti-scratch/kandersen/bin/novocraft_v3/novoalign

# V2.08.03: X86-64 Linux http://www.novocraft.com/main/download.php?filename=V2.08.03/novocraftV2.08.03.gcc.tar.gz
# V2.08.03: Mac          http://www.novocraft.com/main/download.php?filename=V2.08.03/novocraftV2.08.03.MacOSX.tar.gz

# V3.02.02: X86-64 Linux 3.0 Kernel http://www.novocraft.com/main/download.php?filename=V3.02.02/novocraftV3.02.02.Linux3.0.tar.gz
# V3.02.02: X86-64 Linux 2.6 Kernel http://www.novocraft.com/main/download.php?filename=V3.02.02/novocraftV3.02.02.Linux2.6.tar.gz
# V2.02.02: Mac                     http://www.novocraft.com/main/download.php?filename=V3.02.02/novocraftV3.02.02.MacOSX.tar.gz

def get_os_and_version() :
    import os
    uname = os.uname()
    return uname[0], uname[2][:3] # e.g., (Darwin, 13.) or (Linux, 2.6)