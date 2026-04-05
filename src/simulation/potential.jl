function update_electric_field!(∇ϕ, ie_idx, cache, apply_drag)
    # ∇ϕ : mutable array of electric field values (discrete potential gradient over domain)
    #         updated in place. convention in this code: ∇ϕ[i] = -E where E is electric field.
    # ie_idx : index of the intermediate electrode
    # cache : named tuple with cached plasma fields (cell-centered arrays + scalar discharge current)
    #    ji             : ion current density array (A/m^2)
    #    Id             : discharge current scalar (Ref{Float64}, A)
    #    ne             : electron density array (m^-3)
    #    μ              : electron mobility array (m^2/V/s)
    #    ∇pe            : electron pressure gradient array (Pa/m or appropriate unit)
    #    channel_area   : channel cross-section area array (m^2)
    #    νei            : electron-ion collision frequency array (s^-1)
    #    νen            : electron-neutral collision frequency array (s^-1)
    #    νan            : anomalous collision frequency array (s^-1)
    #    avg_ion_vel    : mean ion velocity array (m/s)
    #    avg_neutral_vel: mean neutral velocity array (m/s)
    # apply_drag : Bool; if true include drag corrections from ion/neutrals in E-field update.

    (; ji, Id, Id_L_IE, Id_IE_R, ne, μ, ∇pe, channel_area, νei, νen, νan, avg_ion_vel, avg_neutral_vel) = cache

    if ie_idx > 0
        @printf("  calculating E in double mode\n")
        # calculate E assuming intermediate electrode
        @inbounds for i in eachindex(∇ϕ)
            if i <= ie_idx
                E = ((Id_L_IE[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]
            else
                E = ((Id_IE_R[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]
            end

            if (apply_drag)
                ion_drag = avg_ion_vel[i] * (νei[i] + νan[i]) * me / e
                neutral_drag = avg_neutral_vel[i] * νen[i] * me / e
                E += ion_drag + neutral_drag
            end

            ∇ϕ[i] = -E
        end
    else
        @printf("  calculating E in single mode\n")
        # calcualte E assuming no intermediate electrode (default case)
        @inbounds for i in eachindex(∇ϕ)
            E = ((Id[] / channel_area[i] - ji[i]) / e / μ[i] - ∇pe[i]) / ne[i]

            if (apply_drag)
                ion_drag = avg_ion_vel[i] * (νei[i] + νan[i]) * me / e
                neutral_drag = avg_neutral_vel[i] * νen[i] * me / e
                E += ion_drag + neutral_drag
            end

            ∇ϕ[i] = -E
        end
    end

    return ∇ϕ
end

function integrate_potential!(ϕ, ∇ϕ, grid, V_L)
    # We need to make sure the potential is integrated from the left edge to the right edge,
    # rather than from the left ghost cell center to the right ghost cell center

    # Temorarily replace left and right electric field and grid values with edge values
    EL, ER = ∇ϕ[1], ∇ϕ[end]
    zL, zR = grid.cell_centers[1], grid.cell_centers[end]

    grid.cell_centers[1] = grid.edges[1]
    grid.cell_centers[end] = grid.edges[end]

    ∇ϕ[1] = 0.5 * (EL + ∇ϕ[2])
    ∇ϕ[end] = 0.5 * (ER + ∇ϕ[end - 1])

    # Integrate potential from left to right edge
    cumtrapz!(ϕ, grid.cell_centers, ∇ϕ, V_L)

    # Extrapolate potential to ghost cells
    ϕ[1] = ϕ[1] + (ϕ[1] - ϕ[2])
    ϕ[end] = ϕ[end] + (ϕ[end] - ϕ[end - 1])

    # Replace electric field and cell center values
    grid.cell_centers[1], grid.cell_centers[end] = zL, zR
    ∇ϕ[1], ∇ϕ[end] = EL, ER

    # debug printout
    @printf("  EL: %.3f V/m, ER: %.3f V/m, zL: %.3f cm, zR: %.3f cm\n", EL, ER, zL*100, zR*100)

    return
end


function anode_sheath_potential(params)
    if params.landmark
        return 0.0
    end
    (; anode_bc, cache) = params
    (; ne, ji, channel_area, Tev, Id) = cache

    # Compute anode sheath potential
    @inbounds if anode_bc == :sheath

        Te_sheath_edge = 0.5 * (Tev[1] + Tev[2])
        ne_sheath_edge = 0.5 * (ne[1] + ne[2])
        ce = sqrt(8 * e * Te_sheath_edge / π / me)
        je_sheath = e * ne_sheath_edge * ce / 4

        # discharge current density
        jd = Id[] / channel_area[1]

        # current densities at sheath edge
        ji_sheath_edge = 0.5 * (ji[1] + ji[2])
        je_sheath_edge = jd - ji_sheath_edge

        current_ratio = je_sheath_edge / je_sheath
        if current_ratio ≤ 0.0
            Vs = 0.0
        else
            Vs = -Te_sheath_edge * log(min(1.0, je_sheath_edge / je_sheath))
        end
    else
        Vs = 0.0
    end

    return Vs
end
