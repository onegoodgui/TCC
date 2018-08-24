    program Programa_Trelica

    use Estrutura_Trelica
    use LLS
    use LLSt
    
    implicit none
    ! Variaveis -------------------------------------------------
    integer :: i, iii
    integer, parameter :: uni_dad_estrut = 1080, uni_entrada = 300
    
    !Variaveis LLS
    integer :: n_cant
    type(LLS_var), allocatable :: cant(:)
    
    !Variaveis barra
    integer :: n_barras

    !Corpo do programa principal --------------------------------
    
   ! open (unit=uni_dad_estrut , file="dad_entrada_trelica.txt" , action='read', iostat=iii)   
   !     read(uni_dad_estrut, nml = dad_entrada_trelica)
   ! close (uni_dad_estrut)


    
    open (unit = uni_entrada, file = 'lista_perfis_aco.txt', action='read')
        call armazenar_propriedades_cantoneiras_abas_iguais(n_cant, cant)
    close (unit = uni_entrada)
    
    do i = 1, n_cant
     call LLS_propriedades_geometricas (cant(i))
    end do
        
    call node_barras(n_div, h1, coord)
    call barras_trelica (n_div, h1, coord, n_barras, barra)
    call secao_barras(n_cant, cant, n_barras, barra)
    
    end program Programa_Trelica