#!/bin/bash

vaspT () 
{ 
    if [ ! -f 'OSZICAR' ]; then
        return 1;
    fi;
    if [ -f INCAR ]; then
        local potim=$(grep POTIM INCAR | awk '{print $3}');
    else
        local potim=1;
    fi;
    echo "# Time step: ${potim}" > vasp_md.dat;
    echo "# Step  Temperature Total_energy E_pot E_kin" >> vasp_md.dat;
    awk 'BEGIN {Tsum = 0; Ns = 0};
        /T=/{
            printf "%5d %8.1f %12.6f %12.6f %12.6f\n", $1, $3, $9 + $11, $9, $11;
            Tsum += $3; Ns += 1
        };
        END {print "#", Ns, Tsum/Ns}' OSZICAR >> vasp_md.dat;
    gnuplot  <<EOF
set terminal png size 800,960
set output 'nvt.png' 

set mxtics
set multiplot layout 3,1 rowsfirst

set ylabel "Temperature (K)" font "Candara,20"
set xtics font "Candara,14"
set ytics font "Candara,14"
set lmargin at screen 0.20
p 'vasp_md.dat' u (\$1 * $potim):2 w l lw 3 lc rgb 'black' t ""

set ylabel "PES (eV)" font "Candara,20"
set xtics font "Candara,14"
set ytics font "Candara,14"
set lmargin at screen 0.20
p 'vasp_md.dat' u (\$1 * $potim):4 w l lw 3 lc rgb 'red' t ""

set lmargin at screen 0.20
set xlabel "Time (fs)" font "Candara,20"
set ylabel "Total Energy (eV)" font "Candara,20"
set xtics font "Candara,14"
set ytics font "Candara,14"
p 'vasp_md.dat' u (\$1 * $potim):3 w l lw 3 lc rgb 'blue' t ""

unset multiplot
pause mouse
EOF

    if which feh &>/dev/null; then
        feh -xdF nvt.png;
    fi
}

vaspT
