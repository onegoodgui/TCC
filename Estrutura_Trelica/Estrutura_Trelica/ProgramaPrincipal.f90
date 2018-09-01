    program Programa_Trelica

    use Estrutura_Trelica
    use LLS
    use LLSt
    use sislin_pbMEF
    use Matriz_Rigidez
    
    implicit none
    ! Variaveis -------------------------------------------------
    integer :: i, iii
    integer, parameter :: uni_dad_estrut = 150, uni_entrada = 300, uni_dad_tercas = 500
    
    !Variaveis LLS
    integer :: n_cant
    type(LLS_var), allocatable :: cant(:)
    
    !Variaveis barra
    integer :: n_barras
    real(8) :: rho = 0.000078d0          ! Massa específica do aço [KN/cm³]
    
  
    
    !Corpo do programa principal --------------------------------
    
    !open (unit = uni_dad_estrut, file = 'entrada_trelica.txt', action='read')        
    !    read (uni_dad_estrut, nml = dad_entrada_trelica)
    !close (unit = uni_dad_estrut)
    
    open (unit = uni_dad_tercas, file = 'dados_tercas.txt', action='read')        
        call lista_tercas(uni_dad_tercas, terca)
    close (unit = uni_dad_tercas)
    
    open (unit = uni_entrada, file = 'lista_perfis_aco.txt', action='read')
        call armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)
    close (unit = uni_entrada)
    
    do i = 1, n_cant
     call LLS_propriedades_geometricas (cant(i))
    end do
        
    call node_barras(n_div, h1, num_nos, coord)
    call barras_trelica (n_div, h1, coord, n_barras, barra)
    call secao_barras(n_cant, cant, n_barras, barra)
    call carga_pp_barras (rho, h1, n_div, num_nos, n_barras, barra, cond_cont)
    call carga_tercas_telhas_barras (h1, dist_trelica, n_div, num_nos, barra, terca, cond_cont)
    
    !Matriz de rigidez
    call vetor_Y_vglc (cond_cont, Y, v_glc)
    
    do i = 1, n_barras
        call encontrar_MDC_pbMEF(nne,barra(i)%conectividades(:),MDC)
    end do
    
    call iniciar_MDSB_pbMEF(num_nos,ngln,MDC,MDSB_trelica,err)
    
    allocate(MatRigid(ngln*num_nos,ngln*num_nos))
    MatRigid = 0.d0
    
    do i =1, n_barras
        call matriz_rigidez_local(barra(i)%s%A, E, barra(i)%comprimento, k_local)
        call matriz_trasformacao(barra(i)%node(1)%x,barra(i)%node(1)%y, barra(i)%node(2)%x,barra(i)%node(2)%y, barra(i)%comprimento, T)
        call matriz_rigidez_global(k_local, T, k)
        call matriz_rigidez_estrutura(barra(i)%conectividades(1), barra(i)%conectividades(2), k, MatRigid)
    end do 
    
    call matriz_rigidez_nos_restritos (nglc ,v_glc, MatRigid)
    call eliminacao_gauss(MatRigid, Y)
    call deslocamentos_matriz_rigidez(nglc, v_glc, MatRigid, Y, deslocamentos)
    call forca_axial_barras(barra, n_barras, deslocamentos, forca_axial)
    
    
    end program Programa_Trelica