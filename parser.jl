include("/home/user/Documents/HessianCalc/FFCalc.jl")

function parser(filename::String, Natoms::Int64)
    if filename[end-3:end] == "hess"
        try
            open(filename, "r") do io
                Nfree = Natoms*3
                hessian = zeros(Float64, Nfree, Nfree)
                re5 = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)"
                re4 = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)"
                re3 = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)"
                re2 = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)"
                re1 = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)"
                flag = false
                i = 1
                j = 1
                for ln in eachline(io)
                    if ln === "\$hessian" 
                        flag = true
                    elseif ln === "\$vibrational_frequencies"
                        break
                    elseif flag === true & occursin(re1, ln) === true
                        if length(rsplit(ln)) === 6
                            hess_str = match(re5, ln)
                            hess_val = rsplit(hess_str.match)
                            hess1 = parse(Float64, hess_val[1]); hess2 = parse(Float64, hess_val[2]);
                            hess3 = parse(Float64, hess_val[3]); hess4 = parse(Float64, hess_val[4]);
                            hess5 = parse(Float64, hess_val[5]);
                            hessian[i, j] = hess1; hessian[i, j + 1] = hess2;
                            hessian[i, j + 2] = hess3; hessian[i, j + 3] = hess4;
                            hessian[i, j + 4] = hess5;

                        elseif length(rsplit(ln)) === 5
                            hess_str = match(re4, ln)
                            hess_val = rsplit(hess_str.match)
                            hess1 = parse(Float64, hess_val[1]); hess2 = parse(Float64, hess_val[2]);
                            hess3 = parse(Float64, hess_val[3]); hess4 = parse(Float64, hess_val[4]);
                            hessian[i, j] = hess1; hessian[i, j + 1] = hess2;
                            hessian[i, j + 2] = hess3; hessian[i, j + 3] = hess4;

                        elseif length(rsplit(ln)) === 4
                            hess_str = match(re3, ln)
                            hess_val = rsplit(hess_str.match)
                            hess1 = parse(Float64, hess_val[1]); hess2 = parse(Float64, hess_val[2]);
                            hess3 = parse(Float64, hess_val[3]);
                            hessian[i, j] = hess1; hessian[i, j + 1] = hess2;
                            hessian[i, j + 2] = hess3;

                        elseif length(rsplit(ln)) === 3
                            hess_str = match(re2, ln)
                            hess_val = rsplit(hess_str.match)
                            hess1 = parse(Float64, hess_val[1]); hess2 = parse(Float64, hess_val[2]);
                            hessian[i, j] = hess1; hessian[i, j + 1] = hess2; 

                        elseif length(rsplit(ln)) === 2
                            hess_str = match(re1, ln)
                            hess_val = rsplit(hess_str.match)
                            hess1 = parse(Float64, hess_val[1]);
                            hessian[i, j] = hess1;

                        else
                            continue
                        end

                        j = j + div(i, Nfree) * 5
                        i = i%(Nfree) + 1
                    else
                        continue
                    end
                end
                return(hessian)
            end
        catch
            @warn "File not founded or wrong parsing"
        end

    elseif filename[end-2:end] == "xyz"
        try
            open(filename, "r") do io
                xyz = zeros(Float64, Natoms, 3)
                re = r"([+-]?[0-9]*[.])?[0-9]+\s*[+-]?([0-9]*[.])?[0-9]+\s*[+-]?([0-9]*[.])?[0-9]+(\s*)"
                i = 1
                for ln in eachline(io)
                    if occursin(re, ln) === true
                        m = match(re, ln)
                        coord = rsplit(m.match)
                        xyz[i, 1] = parse(Float64, coord[1])
                        xyz[i, 2] = parse(Float64, coord[2])
                        xyz[i, 3] = parse(Float64, coord[3])
                        i += 1
                    else
                        continue
                    end
                end
                return(xyz)
            end
        catch
            @warn "File not founded or wrong parsing"
        end
    else
        println("Warning!!! Wrong file extension")
    end
end

xyz = parser("/home/user/Documents/HessianCalc/benzene.xyz", 12);
hessian = parser("/home/user/Documents/HessianCalc/benzene.hess", 12);
hessian = round.(hessian; digits=2)
xyz = round.(xyz; digits=2)
system = Hess_data(hessian, xyz)
bonds = [2 1; 3 2; 4 3; 5 4; 6 1; 7 1; 8 2; 9 3; 10 4; 11 5; 12 6; 6 5]
angles  = [1 2 3; 2 3 4;3 4 5;2 1 6;2 1 7;1 2 8;2 3 9;
3 4 10;4 5 11;1 6 12;6 5 11;6 1 7;5 6 12;5 4 10;1 6 5;
4 3 9;4 5 6;3 2 8;]
impropers = [12     6     1     5;
11     5     4     6;
10     4     3     5;
 9     3     2     4;
 8     2     1     3;
 7     1     2     6;]

println("**************************************")
println("*********   BONDS CONSTANTS   ********")
println("**************************************")
for i in 1:size(bonds)[1]
    println(bonds[i, 1], " ", bonds[i, 2], " bond constant:  ", bond_constant(bonds[i, 1] , bonds[i, 2] , system))
end

println("**************************************")
println("*********  ANGLES CONSTANTS   ********")
println("**************************************")
for i in 1:size(angles)[1]
    println(angles[i, 1], " ", angles[i, 2], " ",  angles[i, 3], " angle constant:  ", angle_constant(angles[i, 1] , angles[i, 2], angles[i, 3] , system))
end

println("**************************************")
println("******  IMPROPERS CONSTANTS   ********")
println("**************************************")
for i in 1:size(impropers)[1]
    println(impropers[i, 1], " ", impropers[i, 2], " ",  impropers[i, 3]," ",  impropers[i, 4],  " angle constant:  ", improper_constant(impropers[i, 1], impropers[i, 2], impropers[i, 3], impropers[i, 4] , system))
end
println("**************************************")