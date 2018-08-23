    program Programa_Trelica

    use Estrutura_Trelica
    
    implicit none
    ! Variaveis -------------------------------------------------
    integer :: iii
    integer, parameter :: uni_dad_estrut = 1080

    !Corpo do programa principal --------------------------------
    
   ! open (unit=uni_dad_estrut , file="dad_entrada_trelica.txt" , action='read', iostat=iii)   
   !     read(uni_dad_estrut, nml = dad_entrada_trelica)
   ! close (uni_dad_estrut)

    call node_barras(n_div, h1, coord)
    call barras_trelica (n_div, h1, coord, barra)
    
    end program Programa_Trelica