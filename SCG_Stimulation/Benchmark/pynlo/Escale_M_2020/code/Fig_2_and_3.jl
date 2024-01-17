using LinearAlgebra
using SparseArrays

using ProgressMeter
using Gridap
using Gridap.Geometry
using Gridap.Visualization
using Gridap.ReferenceFEs
using GridapGmsh
using GridapMakie, CairoMakie
using Femwell.Maxwell.Waveguide
using NPZ

function n_LNOI(wavelength::Float64, ray::String ="o")
    if ray == "o"
        if wavelength >= 0.4 && wavelength <=5
            return sqrt(2.6734*wavelength^2/(wavelength^2 - 0.01764) + 1.2290*wavelength^2/(wavelength^2 - 0.05914) + 12.614*wavelength^2/(wavelength^2 - 474.6)+1)
        else
            throw("invalid wavelength for lithium niobate, must be between 0.4-5um")
        end
    elseif ray == "e"
        if wavelength >= 0.4 && wavelength <=5
            return sqrt(2.9804*wavelength^2/(wavelength^2-0.02047)+ 0.5981*wavelength^2/(wavelength^2 -0.0666) + 8.9543 * wavelength^2 / (wavelength^2 - 416.08)+1)
        else
            throw("invalid wavelength for lithium niobate, must be between 0.4-5um")
        end
    end
end

function n_SiO2(wavelength::Float64,type::String ="FusedSilica")
    if type == "FusedSilica"

        if wavelength < 0.21 || wavelength > 6.7
            throw("wavelength provided is $wavelength um, is out of the range for {type}")
        end

        return sqrt( 0.6961663* wavelength^2/(wavelength^2 - 0.0684043^2)+(0.4079426*wavelength^2/(wavelength^2-0.1162414^2))+(0.8974794*wavelength^2/(wavelength^2-9.896161^2))+1)

    elseif type == "flim"
        return 1.45
    end
end


function n_Air(wavelength::Float64)
    return 0.05792105/(238.0185-wavelength^(-2))+0.00167917/(57.362-wavelength^(-2)) + 1
end

function search_by(modes::AbstractArray, type::Function)
    return_mode = modes[1]
    for mode in modes
        if real(type(mode)) >= real(type(return_mode))
            return_mode = mode
        end
    end
    return return_mode
end

function a_eff(mode::Mode)
    e2 =conj(E(mode)) ⋅ E(mode)
    top = sum(∫(e2)Measure(mode))
    bottom = sum(∫(e2*e2)Measure(mode))
    return  top^2/bottom
end

CairoMakie.inline!(true)

wavelength_range = [0.4,1.5] # in um
steps = 70
wavelength_list = LinRange(wavelength_range[1], wavelength_range[2], steps)
neff_list_tm = zeros(0)
neff_list_te = zeros(0)
aeff_list_tm = zeros(0)
aeff_list_te = zeros(0)

# mesh setting
model = GmshDiscreteModel("SCG_Stimulation\\Benchmark\\pynlo\\Escale_M_2020\\code\\mesh.msh")
Ω = Triangulation(model)
labels = get_face_labeling(model)

print("start sweeping")
# Start the sweep
@showprogress dt=1 desc="Computing..." for wavelength in wavelength_list
    # initaliza tensor
    eye = diagonal_tensor(VectorValue([1.0, 1.0, 1.0]))
    LNOI = TensorValue([n_LNOI(wavelength,"e")^2 0 0; 0 n_LNOI(wavelength,"o")^2 0; 0 0 n_LNOI(wavelength,"o")^2])
    epsilons = ["core" => LNOI, "ridge" => LNOI, "buffer" => n_SiO2(wavelength)^2 * eye, "air" => n_Air(wavelength)^2 * eye]

    # Assign tensor
    ε(tag) = Dict(get_tag_from_name(labels, u) => v for (u, v) in epsilons)[tag]
    τ = CellField(get_face_tag(labels, num_cell_dims(model)), Ω)

    # Sort by tm-fraction
    modes = calculate_modes(model, ε ∘ τ, λ = wavelength, num = 2, order = 1)
    mode = search_by(modes, tm_fraction)
    append!(neff_list_tm, n_eff(mode))
    append!(aeff_list_tm, a_eff(mode))

    # Sort by te-fraction
    mode = search_by(modes, te_fraction)
    append!(neff_list_te, n_eff(mode))
    append!(aeff_list_te, a_eff(mode))
    #plot_field(E(modes[1]))
    #plot_mode(mode[1])
end
npzwrite("SCG_Stimulation\\Benchmark\\pynlo\\Escale_M_2020\\code\\data_LNOI_x",  Dict("neff_list_tm" => neff_list_tm, "aeff_list_tm" => aeff_list_tm,"neff_list_te" => neff_list_te, "aeff_list_te" => aeff_list_te, "wls"=>wavelength_list))
# tensor should be in the middle?
#=
CairoMakie.activate!(type = "png")
x = [x * 1e-3 for x = 400:5000]
y = [n_LNOI(y* 1e-3) for y = 400:5000]

f = Figure(resolution = (800, 600));
ax = Axis(f[1, 1], title = "Figure 1", xlabel = "x", ylabel = "y")
lines!(ax, x, y)
sc = display(f)

=#

# use npz module to save into numpy file