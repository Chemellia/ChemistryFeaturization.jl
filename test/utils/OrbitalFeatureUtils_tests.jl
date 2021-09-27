using Test
using SparseArrays: sparsevec
using ..ChemistryFeaturization.Utils.OrbitalFeatureUtils:
    _indexorbital,
    _orbitalindex,
    _orbitalsparse,
    _econf_to_name,
    valenceshell_conf_df

@testset "OrbitalFeatureUtils" begin
    @testset "Index - Orbital" begin
        orbitals =
            SubString.([
                "1s",
                "2s",
                "2p",
                "3s",
                "3p",
                "4s",
                "3d",
                "4p",
                "5s",
                "4d",
                "5p",
                "6s",
                "4f",
                "5d",
                "6p",
                "7s",
                "5f",
                "6d",
                "7p",
                "8s",
                "5g",
                "6f",
                "7d",
                "8p",
                "9s",
            ])
        count = [i for i in Int16.(1:25)]

        @test _indexorbital.(count) == orbitals
        @test _orbitalindex.(orbitals) == count
    end

    @testset "_orbitalsparse" begin
        # in the expected result, the first vector the positions and the second is the values.
        @test _orbitalsparse("1s2") == ([1], [2]) # 1s2
        @test _orbitalsparse("3d5.4s1") == ([7, 6], [5, 1]) # 3d5, 4s1.
    end

    @testset "Element Name - Electronic Configuration" begin

        get_config(el) = _orbitalsparse(getproperty(valenceshell_conf_df[valenceshell_conf_df.Symbol.==el, :][1,:], Symbol("Electronic Structure")))

        function test_econf_to_name(I, V, element)
            @test _econf_to_name(sparsevec(Int16.(I), Int16.(V))) == element
        end

        I_H, V_H = get_config("H")
        @test (I_H, V_H) == ([1], [1])
        test_econf_to_name(I_H, V_H, "H")

        I_He, V_He = get_config("He")
        @test (I_He, V_He) == ([1], [2])
        test_econf_to_name(I_He, V_He, "He")

        I_Li, V_Li = get_config("Li")
        @test (I_Li, V_Li) == ([2], [1])
        test_econf_to_name(I_He, V_He, "He")
    end
end
