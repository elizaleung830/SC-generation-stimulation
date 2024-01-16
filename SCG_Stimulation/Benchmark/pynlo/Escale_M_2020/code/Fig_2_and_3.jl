using LinearAlgebra
using SparseArrays

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using Gridap.ReferenceFEs
using GridapGmsh
using GridapMakie, CairoMakie
using Femwell.Maxwell.Waveguide

function n_LNOI(wavelength::Float32, ray::String ="o")
    if ray == "o"
        if wavelength >= 0.4 && wavelength <=5
            return math.sqrt(2.6734*wavelength^2/(wavelength^2 - 0.01764) + 1.2290*wavelength^2/(wavelength^2 - 0.05914) + 12.614*wavelength^2/(wavelength^2 - 474.6)+1)
        else
            throw("invalid wavelength for lithium niobate, must be between 0.4-5um")
        end
    elseif ray == "e"
        if wavelength >= 0.4 && wavelength <=5
            return math.sqrt(2.9804*wavelength^2/(wavelength^2-0.02047)+ 0.5981*wavelength^2/(wavelength^2 -0.0666) + 8.9543 * wavelength^2 / (wavelength^2 - 416.08)+1)
        else
            throw("invalid wavelength for lithium niobate, must be between 0.4-5um")
        end
    end
end

function n_SiO2(wavelength::Float32,type::String ="FusedSilica")
    if type == "FusedSilica"

        if wavelength < 0.21 || wavelength > 6.7
            throw("wavelength provided is $wavelength um, is out of the range for {type}")
        end

        return np.sqrt( 0.6961663* wavelength^2/(wavelength^2 - 0.0684043^2)+(0.4079426*wavelength^2/(wavelength^2-0.1162414^2))+(0.8974794*wavelength^2/(wavelength^2-9.896161^2))+1)

    elseif type == "flim"
        return 1.45
    end
end


function n_Air(wavelength::Float32)
    return 0.05792105/(238.0185-wavvelength^(-2))+0.00167917/(57.362-wavvelength^(-2)) + 1
end


CairoMakie.inline!(true)
model = GmshDiscreteModel("mesh.msh")
Ω = Triangulation(model)

labels = get_face_labeling(model)

eye = diagonal_tensor(VectorValue([1.0, 1.0, 1.0]))
n2 = 2.302
Δ = 0.005
yig = TensorValue([n2^2 1im*Δ 0; -1im*Δ n2^2 0; 0 0 n2^2])
epsilons = ["core" => yig, "ridge" => 1.95^2 * eye, "buffer" => 1.0^2 * eye, "air" => 1.0^2 * eye]
ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]

τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

modes = calculate_modes(model, ε ∘ τ, λ = 1.3, num = 2, order = 1)

plot_mode(modes[1])
plot_mode(modes[2])
modes

# use npz module to save into numpy file