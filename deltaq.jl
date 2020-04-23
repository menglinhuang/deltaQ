using LinearAlgebra

const mass1 = 1.660539e-27
const mass2 = 9.109383e-31
const bohr = 0.5291772083
const amu_conv = sqrt(1/1822)*0.5291772083

mass_dict = Dict("H"=>1.007, "He"=>4.002, "Li"=>6.941, "Be"=>9.012, "B"=>10.811, "C"=>12.017, "N"=>14.006, "O"=>15.999, "F"=>18.998, "Ne"=>20.179,
"Na"=>22.989, "Mg"=>24.305, "Al"=>26.981, "Si"=>28.085, "P"=>30.973, "S"=>32.065, "Cl"=>35.453, "Ar"=>39.948, "K"=>39.098,
"Ca"=>40.078, "Sc"=>44.955, "Ti"=>47.867, "V"=>50.941, "Cr"=>51.996, "Mn"=>54.938, "Fe"=>55.845, "Co"=>58.933, "Ni"=>58.693,
"Cu"=>63.546, "Zn"=>65.409, "Ga"=>69.723, "Ge"=>72.64, "As"=>74.921, "Se"=>78.96, "Br"=> 79.904, "Kr"=>83.798, "Rb"=>85.467,
"Sr"=>87.62, "Y"=>88.905, "Zr"=>91.224, "Nb"=>92.906, "Mo"=>95.94, "Tc"=>97.9072, "Ru"=>101.07, "Rh"=>102.905, "Pd"=>106.42,
"Ag"=>107.868, "Cd"=>112.411, "In"=>114.818, "Sn"=>118.710, "Sb"=>121.760, "Te"=>127.60, "I"=>126.904, "Xe"=>131.293,
"Cs"=>132.905, "Ba"=>137.327, "La"=>138.905, "Ce"=>140.116, "Pr"=>140.907, "Nd"=>144.242, "Eu"=>151.964, "Hf"=>178.49,
"W"=>183.84, "Pt"=>195.084, "Au"=>196.966, "Hg"=>200.59, "Tl"=>204.383, "Pb"=>207.2)

function change_position(pos_array::Array)
    for value in pos_array
        if value < -0.5
            value = value + 1
        elseif value > 0.5
            value = value - 1
        end
    end
    return pos_array
end

function sum_array(x_array::Array)
    n = 0.0
    for x in x_array
        n += x^2
    end
    return n
end

# read POSCAR-i
file_i = open("POSCAR-i","r")
data_i = readlines(file_i)
close(file_i)
lattice = [parse.(Float64, split(data_i[3])); parse.(Float64, split(data_i[4])); parse.(Float64, split(data_i[5]))]
lattice = permutedims(reshape(lattice,3,3))
atom_type = split(data_i[6])
num_atom_type = length(atom_type)
atom_num = parse.(Int, split(data_i[7]))
tot_atom_num = sum(atom_num)
position_i = []
for i = 1:tot_atom_num
    line = parse.(Float64, split(data_i[i+8]))
    push!(position_i, line)
end

# read POSCAR-f
file_f = open("POSCAR-f","r")
data_f = readlines(file_f)
close(file_f)
position_f = []
for i = 1:tot_atom_num
    line = parse.(Float64, split(data_f[i+8]))
    push!(position_f, line)
end

# define atomic mass
mass = []
for i = 1:num_atom_type
    for j = 1:atom_num[i]
        push!(mass, mass_dict[atom_type[i]])
    end
end

# calculate dR in Cartesian
total_dr = []
sum_dr = 0.0
for i = 1:tot_atom_num
    dr = permutedims(change_position(position_f[i]-position_i[i]))*lattice
    sqrt(dr[1]^2+dr[2]^2+dr[3]^2) > 3 && throw(ErrorException("Too large dR!"))
    push!(total_dr, dr)
    global sum_dr += dr[1]^2+dr[2]^2+dr[3]^2
end

vector = total_dr/sqrt(sum_dr)  # define eigenvector

dq = 0.0
norm = 0.0
for i = 1:tot_atom_num
    global dq += sqrt(mass[i]*mass1/mass2)*dot(vector[i],total_dr[i]/bohr)
    global norm += sum_array(vector[i])
end
dq_amu = dq*amu_conv
println("Check normalization: $norm")
println("delta_Q = $dq_amu amu^1/2*Angs")
