using LinearAlgebra

# Calculation of Intramolecular Force Fields from Second-Derivative Tensors
# Jorge M. Seminario

# Data struct 
# Orca uses energy in Hartree units and length in Bohr
# 1 Hartee = 2625.5 kJ/mol; 1 Bohr = 0.0529177249 nm
# 1 [Hartee/Bohr^2] = 937582.935092 [kJ/mol/nm^2]
# 1 [Hartee/Bohr^2] = 2240.8815787 [kcal/mol/A^2]
#  xyz in angstrom; 1 A = 0.1 nm
struct Hess_data
    Hess::Matrix{Float64}
    xyz::Matrix{Float64}
    function Hess_data(Hess::Matrix{Float64}, xyz::Matrix{Float64})
        new(Hess*937582.935092 , 0.1*xyz)
    end
end

function bond_constant(atom1::Int64, atom2::Int64, MolSystem::Hess_data)
    # getting atom1-atom2 bond vector
    bond_vector = MolSystem.xyz[atom2, :] - MolSystem.xyz[atom1, :]
    bond_vector = bond_vector/norm(bond_vector)

    # getting submatrice of atom1-atom2 interaction
    HessAB = MolSystem.Hess[(3*atom1-2):(3*atom1), (3*atom2-2):3*atom2]
    eig, vec = eigen(-HessAB)
    
    # calculating bond force constant
    k = eig[1]*abs(dot(bond_vector, vec[:,1])) + eig[2]*abs(dot(bond_vector, vec[:,2])) + eig[3]*abs(dot(bond_vector, vec[:,3]))


    return(round(real(k), digits=3))

end

function angle_constant(atom1::Int64, atom2::Int64, atom3::Int64, MolSystem::Hess_data)
    # getting atom1-atom2; atom3-atom2 bond unit vectors ant their norms
    bond12_vector = MolSystem.xyz[atom2, :] - MolSystem.xyz[atom1, :]
    bond32_vector = MolSystem.xyz[atom2, :] - MolSystem.xyz[atom3, :]
    bond12_dist = norm(bond12_vector)
    bond32_dist = norm(bond32_vector)
    bond12_vector = bond12_vector/bond12_dist
    bond32_vector = bond32_vector/bond32_dist

    # getting submatrices of atom1-atom2 and atom2-atom3 interactions
    # calculating their eigen -values and -vectors
    Hess12 = MolSystem.Hess[(3*atom1-2):(3*atom1), (3*atom2-2):3*atom2]
    Hess32 = MolSystem.Hess[(3*atom3-2):(3*atom3), (3*atom2-2):3*atom2]
    eig12, vec12 = eigen(-Hess12)
    eig32, vec32 = eigen(-Hess32)

    # calculating direction vectors of small angle displacements 
    normN_vector = cross(bond32_vector, bond12_vector)/norm(cross(bond32_vector, bond12_vector))
    disp12_vector = cross(normN_vector, bond12_vector)/norm(cross(normN_vector, bond12_vector))
    disp32_vector = cross(bond32_vector, normN_vector)/norm(cross(bond32_vector, normN_vector))

    # calculating angle force constant
    tempK1 = eig12[1]*abs(dot(disp12_vector, vec12[:,1])) + eig12[2]*abs(dot(disp12_vector, vec12[:,2])) + eig12[3]*abs(dot(disp12_vector, vec12[:,3]))
    tempK2 = eig32[1]*abs(dot(disp32_vector, vec32[:,1])) + eig32[2]*abs(dot(disp32_vector, vec32[:,2])) + eig32[3]*abs(dot(disp32_vector, vec32[:,3]))
    k = 1/(bond12_dist^2*tempK1) + 1/(bond32_dist^2*tempK2)
    k_theta = 1/k/2

    return(round(real(k_theta); digits=3))
end

function dihedral_constant(atom1::Int64, atom2::Int64, atom3::Int64, atom4::Int64,  MolSystem::Hess_data)
    # getting atom2-atom1; atom2-atom3; atom3-atom4; atom3-atom2 bond unit vectors ant their norms
    bond12_vector = MolSystem.xyz[atom2, :] - MolSystem.xyz[atom1, :]
    bond32_vector = MolSystem.xyz[atom2, :] - MolSystem.xyz[atom3, :]
    bond43_vector = MolSystem.xyz[atom3, :] - MolSystem.xyz[atom4, :]
    bond12_dist = norm(bond12_vector)
    bond32_dist = norm(bond32_vector)
    bond43_dist = norm(bond43_vector)
    bond12_vector = bond12_vector/bond12_dist
    bond32_vector = bond32_vector/bond32_dist
    bond43_vector = bond43_vector/bond43_dist
    bond23_vector = -bond32_vector/bond32_dist

    # getting submatrices of atom1-atom2 and atom3-atom4 interactions
    # calculating their eigen -values and -vectors
    Hess12 = MolSystem.Hess[(3*atom1-2):(3*atom1), (3*atom2-2):3*atom2]
    Hess34 = MolSystem.Hess[(3*atom3-2):(3*atom3), (3*atom4-2):3*atom4]
    eig12, vec12 = eigen(-Hess12)
    eig34, vec34 = eigen(-Hess34)

    # calculating unit normal vectors for planes atom1-atom2-atom3(ABC) and atom2-atom3-atom4(BCD)
    normNABC_vector = cross(bond32_vector, bond12_vector)
    normNBCD_vector = cross(bond43_vector, bond23_vector)
    normNABC_dist = norm(normNABC_vector)
    normNBCD_dist = norm(normNBCD_vector)
    normNABC_vector = normNABC_vector/normNABC_dist
    normNBCD_vector = normNBCD_vector/normNBCD_dist

    # calculating dihedral angle force constant
    tempK1 = eig12[1]*abs(dot(normNABC_vector, vec12[:,1])) + eig12[2]*abs(dot(normNABC_vector, vec12[:,2])) + eig12[3]*abs(dot(normNABC_vector, vec12[:,3]))
    tempK2 = eig34[1]*abs(dot(normNBCD_vector, vec34[:,1])) + eig34[2]*abs(dot(normNBCD_vector, vec34[:,2])) + eig34[3]*abs(dot(normNBCD_vector, vec34[:,3]))
    k = 1/((bond12_dist*normNABC_dist)^2*tempK1) + 1/((bond34_dist*normNBCD_dist)^2*tempK2)
    k_phi = 1/k

    return(abs(real(k_phi)))
end

function improper_constant(atom1::Int64, atom2::Int64, atom3::Int64, atom4::Int64,  MolSystem::Hess_data)
    # getting atom1-atom2; atom1-atom3; atom1-atom4 and atom4-atom3; atom2-atom3 bond unit vectors ant their norms
    # getting coordinates of central atom1 and atom2
    Atom1_point = MolSystem.xyz[atom1, :]
    Atom2_point = MolSystem.xyz[atom2, :]
    bond12_vector = Atom2_point - Atom1_point
    bond13_vector = MolSystem.xyz[atom3, :] - Atom1_point
    bond14_vector = MolSystem.xyz[atom4, :] - Atom1_point
    bond43_vector = MolSystem.xyz[atom3, :] - MolSystem.xyz[atom4, :]
    bond23_vector = MolSystem.xyz[atom3, :] - Atom2_point
    bond12_dist = norm(bond12_vector)
    bond13_dist = norm(bond13_vector)
    bond14_dist = norm(bond14_vector)
    bond43_dist = norm(bond43_vector)
    bond23_dist = norm(bond23_vector)
    bond12_vector = bond12_vector/bond12_dist
    bond13_vector = bond13_vector/bond13_dist
    bond14_vector = bond14_vector/bond14_dist
    bond43_vector = bond43_vector/bond43_dist
    bond23_vector = bond23_vector/bond23_dist


    # getting submatrices of atom1-atom2 and atom1-atom3 and atom1-atom4 interactions
    # calculating their eigen -values and -vectors
    Hess12 = MolSystem.Hess[(3*atom1-2):(3*atom1), (3*atom2-2):3*atom2]
    Hess13 = MolSystem.Hess[(3*atom1-2):(3*atom1), (3*atom3-2):3*atom3]
    Hess14 = MolSystem.Hess[(3*atom1-2):(3*atom1), (3*atom4-2):3*atom4]
    eig12, vec12 = eigen(-Hess12)
    eig13, vec13 = eigen(-Hess13)
    eig14, vec14 = eigen(-Hess14)

    # calculating unit normal vectors for planes atom2-atom3-atom4(BCD)
    normNBCD_vector = cross(bond43_vector, bond23_vector)
    normNBCD_dist = norm(normNBCD_vector)
    normNBCD_vector = normNBCD_vector/normNBCD_dist

    # calculating improper dihedral angle force constant 

    # calculating in
    tempK1 = eig12[1]*abs(dot(normNBCD_vector , vec12[:,1])) + eig12[2]*abs(dot(normNBCD_vector , vec12[:,2])) + eig12[3]*abs(dot(normNBCD_vector , vec12[:,3]))
    tempK2 = eig13[1]*abs(dot(normNBCD_vector , vec13[:,1]))+ eig13[2]*abs(dot(normNBCD_vector , vec13[:,2])) + eig13[3]*abs(dot(normNBCD_vector , vec13[:,3]))
    tempK3 = eig14[1]*abs(dot(normNBCD_vector , vec14[:,1])) + eig14[2]*abs(dot(normNBCD_vector , vec14[:,2])) + eig14[3]*abs(dot(normNBCD_vector , vec14[:,3]))

    # calculating projection of altitude of triangle ABC into plane BCD
    # determinating altitude vector of altitude through system of linear equations
    # Statement 1: equation of line BC; Statement 2: dot product is 0 for (BC*altitude),  MatrCoef*r = x
    # BC is bond23_vector; Atom2_point is B point; Atom1_point is A point
     
    MatrCoef = [1/bond23_vector[1]   -1/bond23_vector[2]           0;
                          0           1/bond23_vector[2]   -1/bond23_vector[3];
                bond23_vector[1]        bond23_vector[2]    bond23_vector[3]]

    x = [(Atom2_point[1]/bond23_vector[1] - Atom2_point[2]/bond23_vector[2]);
         (Atom2_point[2]/bond23_vector[2] - Atom2_point[3]/bond23_vector[3]);
        -(dot(Atom1_point, bond23_vector))]


    HABC_vector = Atom1_point - MatrCoef\x
    normHABCD_dist = norm(HABC_vector - HABC_vector*dot(HABC_vector, normNBCD_vector)) 
    k_omega = normHABCD_dist^2*(tempK1 + tempK2 + tempK3)

    return(abs(real(k_omega)))
end