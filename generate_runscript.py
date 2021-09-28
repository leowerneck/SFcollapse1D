import os,sys

def run_string(eta,optimization_string):
    return optimization_string+"./SFcollapse1D 320 16 5.3162 0.08 "+eta

def generate_runscript():
    if sys.platform == "linux" or sys.platform == "linux2":

        # Use taskset to improve code's performance
        import multiprocess as mp

        # Get only physical cores
        N_cores = int(mp.cpu_count()/2)

        # Write string
        optimization_string = "taskset -c 0"
        for i in range(1,N_cores):
            optimization_string += ","+str(i)
        optimization_string += " "

    elif sys.platform == "darwin":

        # Mac OS does not support taskset. Do nothing.
        optimization_string = ""

    else:
        # Windows
        print("Windows detected. runscript.sh will not be generated. Please run the code manually.")
        return

    # Set eta_weak and eta_strong
    eta_weak   = "0.3364266156435"
    eta_strong = "0.3364266156436"

    # Initialize file with the bash environment
    filestring  = "#!/bin/bash\n"
    filestring += run_string(eta_weak,optimization_string)+"\n"
    filestring += "mv out_central_values.dat out_weak.dat\n"
    filestring += run_string(eta_strong,optimization_string)+"\n"
    filestring += "mv out_central_values.dat out_strong.dat\n"
    filestring += "python generate_plot.py"

    with open("runscript.sh","w") as file:
        # Write string to file
        file.write(filestring)

    os.system("chmod 755 runscript.sh")

if __name__ == '__main__':
    generate_runscript()
