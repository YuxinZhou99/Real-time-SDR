# example.gnuplut : configuration for plotting (change as needed)

reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,1.0  # set two plots for this figure

set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2

# freq domain (Fourier)
# set ylabel 'Imaginary'              # set y-axis label
# set xlabel 'Real'               # set x-axis label
# set yrange [-1.2:1.2]                    # set y plot range
# set xrange [300:400]                       # set x plot range
# plot '../data/rrc_output.dat' with line
# plot '../data/rrc_output.dat' with lines linestyle 1, \
 '../data/rrc_output1.dat' with lines linestyle 2

# freq domain (Fourier)
set ylabel 'Imaginary'              # set y-axis label
set xlabel 'Real'               # set x-axis label
set yrange [-0.6:0.7]                    # set y plot range
set xrange [0:610]                       # set x plot range
plot '../data/rrc_output1.dat' with line
# plot '../data/rrc_output1.dat' with lines linestyle 1, \
 '../data/test.dat' with lines linestyle 2

# scatter
set ylabel 'Imaginary'              # set y-axis label
set xlabel 'Real'               # set x-axis label
set yrange [-0.6:0.6]                    # set y plot range
set xrange [-1:1]                       # set x plot range
plot '../data/rrc_output.dat' with points pt 7 ps 1


unset multiplot
