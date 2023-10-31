module macmods_mod

    !-----------------------------------------------------------------------
    !  The following are used extensively in this ecosys, so are used at
    !  the module level. The use statements for variables that are only needed
    !  locally are located at the module subprogram level.
    !-----------------------------------------------------------------------
    use marbl_kinds_mod, only : int_kind
    use marbl_kinds_mod, only : r8

    use marbl_constants_mod, only : c0
    use marbl_constants_mod, only : c1
    use marbl_constants_mod, only : c2
    use marbl_constants_mod, only : c3
    use marbl_constants_mod, only : c4
    use marbl_constants_mod, only : c10
    use marbl_constants_mod, only : c1000
    use marbl_constants_mod, only : p001
    use marbl_constants_mod, only : p5
    use marbl_constants_mod, only : pi


    subroutine MACMODS_init()

    end subroutine MACMODS_init

    ! MACMODS tracers, we don't move!
    ['B', 'QN', 'QP', 'QFe', 'Gave', 'Dave']  ! variables that are carried as tracers (in/out during model calculation)

    output_vars = ['d_B', 'd_QN', 'd_QP', 'd_QFe', 'Growth2', 'd_Be',
                   'd_Bm', 'd_Ns', 'harv', 'GRate', 'B_N', 'n_harv',
                   'min_lim', 'gQ', 'gT', 'gE', 'gH']


    subroutine macmods_load_parameters()
    end subroutine macmods_load_parameters

    subroutine macmods_init_timing()
    end subroutine macmods_init_timing

    ! register tracers, or allocate if local
    ! register outputs, or allocate if local
    ! register internal mem? Seeding/harvest tables?
    subroutine macmods_init_memory()
    end subroutine

    ! init surface and other forcing data structures
    subroutine macmods_init_forcing()
    end subroutine

    ! use run parameters to construct seeding and harvest setup
    subroutine macmods_init_seeding_harvest()
    end subroutine

    ! all the stuff in compute currently, plus:
    ! TODO: establish NH4 uptake and params
    ! TODO: establish PO4 uptake and params
    ! TODO: establish Fe uptake and params
    ! TODO: establish particulate export - what makes it into marbl POC and which goes straight to bottom cell or goes poof
    ! TODO: establish DOM export, stoichiometry of export
    ! TODO: establish stoichiometry parameters/ratios and machinery for limiting growth
    ! TODO: partition death into a macro-grazing, respiration/background death, breakage/POC,

    subroutine macmods_tendency_compute()
    end subroutine

    ! TODO: how to write to outputs
    subroutine macmods_report()
    end subroutine macmods_report

    ! ecosys tracers,


    ! need have


    ! outputs: tendencies for tracers:


end module macmods_mod