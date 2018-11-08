module Estrutura_Trelica

    use LLSt
    use LLScross
    ! módulo para calculo dos esforços numa estrutura de treliça trapezoidal/triangular, de número de
    implicit none

    
    ! Variáveis TYPE -------------------------------------------------------------------------------------
    type coordenadas
        real(8) :: x   ! [cm]
        real(8) :: y   ! [cm]
    end type
    
    type node
        type(coordenadas), allocatable :: no(:)
    end type
    
    type cond_contorno       
        real(8) :: carga_pp(2)              ! carga (1) -> carga na direção x / carga(2) -> carga na direção y
        real(8) :: carga_sobrecarga(2)      ! carga (1) -> carga na direção x / carga(2) -> carga na direção y
        real(8) :: carga_telha_telhado(2)   ! carga (1) -> carga na direção x / carga(2) -> carga na direção y 
        real(8) :: carga_vento(2)           ! carga (1) -> carga na direção x / carga(2) -> carga na direção y 
        real(8) :: glc(2)                   ! glc(1) -> grau de liberdade da translação na direção do eixo x / glc(2) -> grau de liberdade da translação na direção do eixo y
    end type
    
    type terca_list
        character(20) :: nome
        real(8) :: carga_linear        ! [KN/cm]
        real(8) :: A                   ! [cm²]
    end type
    

    type barra_trelica
        type(coordenadas) :: node(2)
        integer :: conectividades(2)
        real(8) :: comprimento        ! [cm]
        type(LLSt_var) :: s
        real(8) :: peso
        character(7) :: tipo
        real(8) :: esbeltez
    end type
        
    
    type comb_acoes
        character(30) :: nome
        type(coordenadas), allocatable :: carga_nodal(:)
    end type
    
    ! Variaveis para montagem da estrutura ---------------------------------------------------------------

    integer :: theta                                           ! numero total de nós da estrutura
    integer :: num_nos
    type(node) :: coord
    real(8) :: L                                               ! vão horizontal entre as extremidades do pórtico [cm]
    real(8) :: h1                                              ! altura do montante de extremidade [cm]
    real(8) :: dist_trelica                                    ! distância do vão entre dois pórticos consecutivos [cm]
    integer :: n_div                                           ! divisões do comprimento referente a L/2 [cm]
    character(25) :: inclinacao_diagonais                      !
    type(terca_list), allocatable :: terca(:)

    real(8), parameter :: pi = 4.d0*atan(1.d0)
    
    contains
    
    
    
    !***********************************************************************************************************************************************
    subroutine node_barras(n_div, h1, num_nos, coord)  !Subrotina que numera os nós da treliça e calcula as coordenadas de cada nó na sua respectiva numeração
    !***********************************************************************************************************************************************
    integer(4), intent(in) :: n_div           ! número de divisões de L/2
    real(8), intent(in) :: h1                 ! altura do montante de extremidade da cobertura treliçada
    integer, intent(out) :: num_nos           ! numero total de nós da estrutura
    type(node), intent(out) :: coord          ! variável de armazenamento das coordenadas dos nós
    
    integer :: i = 0, n = 0
    
    
        if (h1>0) then
            
            !---------COBERTURA TRAPEZOIDAL -----------
            num_nos = 4*n_div + 2
            allocate(coord%no(num_nos))
             
                ! atribuição de valores para as coordenadas x dos nós
            do i = 1, 2*n_div+1
                coord%no(2*(i)-1)%x  = (i-1)*L/(2*n_div)
                coord%no(2*(i))%x = (i-1)*L/(2*n_div)
            end do
                
                ! atribuição de valores para as coordenadas y dos nós
             do i = 1, 2*n_div+1
                 
                 if(2*(i)-1 > 2*n_div+1) then
                     n=n+1
                         coord%no(2*(i)-1)%y = 0
                         coord%no(2*(i))%y = coord%no(2*(i) - n*4)%y
                         cycle
                 end if
                 
                coord%no(2*(i)-1)%y = 0.d0 
                coord%no(2*(i))%y = h1 + (i-1)*L/(2*n_div)*tan(theta*pi/180)
             end do
             n=0
        else
            !---------COBERTURA TRIANGULAR -----------
            num_nos = 4*n_div
            allocate(coord%no(num_nos))
            
                ! coordenadas do primeiro e último nó
            coord%no(1)%x  = 0
            coord%no(4*n_div)%x = L
            coord%no(1)%y  = 0
            coord%no(4*n_div)%y = 0
            
                ! atribuição de valores para as coordenadas x dos nós
            do i = 1, 2*n_div - 1
                coord%no(2*(i))%x  = (i)*L/(2*n_div)
                coord%no(2*(i)+1)%x = (i)*L/(2*n_div)
            end do
            
                ! atribuição de valores para as coordenadas y dos nós
             do i = 1, 2*n_div-1
                 if(2*(i) > 2*n_div) then
                     n=n+1
                         coord%no(2*(i))%y = 0
                         coord%no(2*(i)+1)%y = coord%no(2*(i)+1 - n*4)%y
                         cycle
                 end if
                 
                coord%no(2*(i))%y = 0.d0 
                coord%no(2*(i)+1)%y = h1 + L/(2*n_div)*(i)*tan(theta*pi/180)
             end do
        end if
        
                
    end subroutine
    
    !***************************************************************************************************************************************
    subroutine barras_trelica (n_div, h1, coord, inclinacao_diagonais, n_barras, barra)
    !***************************************************************************************************************************************
        integer(4), intent(in) :: n_div                                ! nº de divisões de L/2
        real(8), intent(in) :: h1                                      ! altura do montante de extremidade da cobertura
        type(node), intent(in) :: coord                                ! vetor com os nós e suas respectivas coordenadas
        character(25), intent(in) :: inclinacao_diagonais
        integer, intent(out) :: n_barras                               ! número de barras da estrutura treliçada
        type(barra_trelica), intent(out), allocatable :: barra(:)      ! vetor com as barras, suas conectividades e seu comprimento
        
        
        integer :: n=0, i=0, k=0                                      ! contadores
        
        if(h1>0) then
            !---------COBERTURA TRAPEZOIDAL -----------
            n_barras = 8*n_div+1
            allocate(barra(n_barras))
            
            ! barras à esquerda em relação ao eixo de simetria ------------------            
            do i = 1, 2*n_div
                n = i+1

                                
                barra(2*i-1)%conectividades(1) = i
                barra(2*i-1)%conectividades(2) = i+1
                barra(2*i)%conectividades(1) = i
                barra(2*i)%conectividades(2) = i+2

                
                barra(2*i-1)%node(1)%x = coord%no(i)%x
                barra(2*i-1)%node(1)%y = coord%no(i)%y
                barra(2*i-1)%node(2)%x = coord%no(i+1)%x
                barra(2*i-1)%node(2)%y = coord%no(i+1)%y
            
                barra(2*i)%node(1)%x = coord%no(i)%x
                barra(2*i)%node(1)%y = coord%no(i)%y
                barra(2*i)%node(2)%x = coord%no(i+2)%x
                barra(2*i)%node(2)%y = coord%no(i+2)%y
            !-------------------------------------------------------------    
                
            end do
            !barra do meio --------------------------------------
            
            barra(4*n_div+1)%conectividades(1) = i
            barra(4*n_div+1)%conectividades(2) = i+1
            
            barra(4*n_div+1)%node(1)%x = coord%no(i)%x
            barra(4*n_div+1)%node(1)%y = coord%no(i)%y
            
            barra(4*n_div+1)%node(2)%x = coord%no(i+1)%x
            barra(4*n_div+1)%node(2)%y = coord%no(i+1)%y
            !---------------------------------------------------
            
            i=i+2
            k = 0
            
            ! barras da direita em relação ao eixo de simetria --------------------------------
            do i = i, 4*n_div + 2
                
                if (mod(i,2) == 0) then
                
                barra(4*(n_div+k)+3)%conectividades(1) = i-3
                barra(4*(n_div+k)+3)%conectividades(2) = i
                barra(4*(n_div+k)+4)%conectividades(1) = i-2
                barra(4*(n_div+k)+4)%conectividades(2) = i
                
                barra(4*(n_div+k)+3)%node(1)%x = coord%no(i-3)%x
                barra(4*(n_div+k)+3)%node(1)%y = coord%no(i-3)%y
                barra(4*(n_div+k)+3)%node(2)%x = coord%no(i)%x
                barra(4*(n_div+k)+3)%node(2)%y = coord%no(i)%y
            
                barra(4*(n_div+k)+4)%node(1)%x = coord%no(i-2)%x
                barra(4*(n_div+k)+4)%node(1)%y = coord%no(i-2)%y
                barra(4*(n_div+k)+4)%node(2)%x = coord%no(i)%x
                barra(4*(n_div+k)+4)%node(2)%y = coord%no(i)%y
                
                k=k+1
                cycle
                end if
            
                
                barra(4*(n_div+k)+2)%conectividades(1) = i-2
                barra(4*(n_div+k)+2)%conectividades(2) = i
                barra(4*(n_div+k+1)+1)%conectividades(1) = i
                barra(4*(n_div+k+1)+1)%conectividades(2) = i+1
                
                barra(4*(n_div+k)+2)%node(1)%x = coord%no(i-2)%x
                barra(4*(n_div+k)+2)%node(1)%y = coord%no(i-2)%y
                barra(4*(n_div+k)+2)%node(2)%x = coord%no(i)%x
                barra(4*(n_div+k)+2)%node(2)%y = coord%no(i)%y
            
                barra(4*(n_div+k+1)+1)%node(1)%x = coord%no(i)%x
                barra(4*(n_div+k+1)+1)%node(1)%y = coord%no(i)%y
                barra(4*(n_div+k+1)+1)%node(2)%x = coord%no(i+1)%x
                barra(4*(n_div+k+1)+1)%node(2)%y = coord%no(i+1)%y
            
            end do
            
            if(inclinacao_diagonais == "/\") then
                do i = 1, 2*n_div
                    if(3 +4*(i-1) < 4*n_div +1) then
                        barra(3 +4*(i-1))%conectividades(1) = barra(3 +4*(i-1))%conectividades(1) - 1
                        barra(3 +4*(i-1))%conectividades(2) = barra(3 +4*(i-1))%conectividades(2) + 1
                        
                        barra(3 +4*(i-1))%node(1)%x = coord%no(barra(3 +4*(i-1))%conectividades(1))%x
                        barra(3 +4*(i-1))%node(1)%y = coord%no(barra(3 +4*(i-1))%conectividades(1))%y
                        barra(3 +4*(i-1))%node(2)%x = coord%no(barra(3 +4*(i-1))%conectividades(2))%x
                        barra(3 +4*(i-1))%node(2)%y = coord%no(barra(3 +4*(i-1))%conectividades(2))%y
                        
                    else if (3 +4*(i-1) > 4*n_div +1) then
                        barra(3 +4*(i-1))%conectividades(1) = barra(3 +4*(i-1))%conectividades(1) + 1
                        barra(3 +4*(i-1))%conectividades(2) = barra(3 +4*(i-1))%conectividades(2) - 1
                                                                    
                        barra(3 +4*(i-1))%node(1)%x = coord%no(barra(3 +4*(i-1))%conectividades(1))%x
                        barra(3 +4*(i-1))%node(1)%y = coord%no(barra(3 +4*(i-1))%conectividades(1))%y
                        barra(3 +4*(i-1))%node(2)%x = coord%no(barra(3 +4*(i-1))%conectividades(2))%x
                        barra(3 +4*(i-1))%node(2)%y = coord%no(barra(3 +4*(i-1))%conectividades(2))%y
                    end if
                end do
                
                    
                else if (inclinacao_diagonais == "mista") then
                    barra(3)%conectividades(1) = barra(3)%conectividades(1) - 1
                    barra(3)%conectividades(2) = barra(3)%conectividades(2) + 1
                    
                        barra(3)%node(1)%x = coord%no(barra(3)%conectividades(1))%x
                        barra(3)%node(1)%y = coord%no(barra(3)%conectividades(1))%y
                        barra(3)%node(2)%x = coord%no(barra(3)%conectividades(2))%x
                        barra(3)%node(2)%y = coord%no(barra(3)%conectividades(2))%y
                        
                    barra(4*(n_div-1)+3)%conectividades(1) = barra(4*(n_div-1)+3)%conectividades(1) - 1
                    barra(4*(n_div-1)+3)%conectividades(2) = barra(4*(n_div-1)+3)%conectividades(2) + 1
                    
                        barra(4*(n_div-1)+3)%node(1)%x = coord%no(barra(4*(n_div-1)+3)%conectividades(1))%x
                        barra(4*(n_div-1)+3)%node(1)%y = coord%no(barra(4*(n_div-1)+3)%conectividades(1))%y
                        barra(4*(n_div-1)+3)%node(2)%x = coord%no(barra(4*(n_div-1)+3)%conectividades(2))%x
                        barra(4*(n_div-1)+3)%node(2)%y = coord%no(barra(4*(n_div-1)+3)%conectividades(2))%y
                    
                    barra(4*(n_div)+3)%conectividades(1) = barra(4*(n_div)+3)%conectividades(1) + 1
                    barra(4*(n_div)+3)%conectividades(2) = barra(4*(n_div)+3)%conectividades(2) - 1
                    
                        barra(4*(n_div)+3)%node(1)%x = coord%no(barra(4*(n_div)+3)%conectividades(1))%x
                        barra(4*(n_div)+3)%node(1)%y = coord%no(barra(4*(n_div)+3)%conectividades(1))%y
                        barra(4*(n_div)+3)%node(2)%x = coord%no(barra(4*(n_div)+3)%conectividades(2))%x
                        barra(4*(n_div)+3)%node(2)%y = coord%no(barra(4*(n_div)+3)%conectividades(2))%y
                        
                    barra(8*(n_div)-1)%conectividades(1) = barra(8*(n_div)-1)%conectividades(1) + 1
                    barra(8*(n_div)-1)%conectividades(2) = barra(8*(n_div)-1)%conectividades(2) - 1
                    
                        barra(8*(n_div)-1)%node(1)%x = coord%no(barra(8*(n_div)-1)%conectividades(1))%x
                        barra(8*(n_div)-1)%node(1)%y = coord%no(barra(8*(n_div)-1)%conectividades(1))%y
                        barra(8*(n_div)-1)%node(2)%x = coord%no(barra(8*(n_div)-1)%conectividades(2))%x
                        barra(8*(n_div)-1)%node(2)%y = coord%no(barra(8*(n_div)-1)%conectividades(2))%y
                else
            end if
            
                
        else
            !---------COBERTURA TRIANGULAR -----------
            n_barras = 8*n_div-3
            allocate(barra(n_barras))
            
            do i = 1, 2*n_div-1
                                
                barra(2*i-1)%conectividades(1) = i
                barra(2*i-1)%conectividades(2) = i+1
                barra(2*i)%conectividades(1) = i
                barra(2*i)%conectividades(2) = i+2

                
                barra(2*i-1)%node(1)%x = coord%no(i)%x
                barra(2*i-1)%node(1)%y = coord%no(i)%y
                barra(2*i-1)%node(2)%x = coord%no(i+1)%x
                barra(2*i-1)%node(2)%y = coord%no(i+1)%y
            
                barra(2*i)%node(1)%x = coord%no(i)%x
                barra(2*i)%node(1)%y = coord%no(i)%y
                barra(2*i)%node(2)%x = coord%no(i+2)%x
                barra(2*i)%node(2)%y = coord%no(i+2)%y
            end do
            
            !barra do meio --------------------------------------
            
            barra(4*n_div-1)%conectividades(1) = i+1
            barra(4*n_div-1)%conectividades(2) = i+2
            
            barra(4*n_div-1)%node(1)%x = coord%no(i+1)%x
            barra(4*n_div-1)%node(1)%y = coord%no(i+1)%y
            
            barra(4*n_div-1)%node(2)%x = coord%no(i+2)%x
            barra(4*n_div-1)%node(2)%y = coord%no(i+2)%y
            !---------------------------------------------------
            
            i=i+2
            k = 0
            
            ! barras da direita em relação ao eixo de simetria --------------------------------
            do i = i, 4*n_div
                
                if (mod(i,2) == 1) then
                
                    barra(4*(n_div+k)+2)%conectividades(1) = i-2
                    barra(4*(n_div+k)+2)%conectividades(2) = i
                    barra(4*(n_div+k)+1)%conectividades(1) = i-1
                    barra(4*(n_div+k)+1)%conectividades(2) = i
                
                    barra(4*(n_div+k)+2)%node(1)%x = coord%no(i-2)%x
                    barra(4*(n_div+k)+2)%node(1)%y = coord%no(i-2)%y
                    barra(4*(n_div+k)+2)%node(2)%x = coord%no(i)%x
                    barra(4*(n_div+k)+2)%node(2)%y = coord%no(i)%y
            
                    barra(4*(n_div+k)+1)%node(1)%x = coord%no(i-1)%x
                    barra(4*(n_div+k)+1)%node(1)%y = coord%no(i-1)%y
                    barra(4*(n_div+k)+1)%node(2)%x = coord%no(i)%x
                    barra(4*(n_div+k)+1)%node(2)%y = coord%no(i)%y
                
                    k=k+1
                
                    cycle
                
                end if
                
                if (i == 4*n_div) then
                    barra(8*(n_div)-4)%conectividades(1) = 4*n_div-2
                    barra(8*(n_div)-4)%conectividades(2) = 4*n_div
                    barra(8*(n_div)-3)%conectividades(1) = 4*n_div-1
                    barra(8*(n_div)-3)%conectividades(2) = 4*n_div
                    
                    barra(8*(n_div)-4)%node(1)%x = coord%no(4*n_div-2)%x
                    barra(8*(n_div)-4)%node(1)%y = coord%no(4*n_div-2)%y
                    barra(8*(n_div)-4)%node(2)%x = coord%no(4*n_div)%x
                    barra(8*(n_div)-4)%node(2)%y = coord%no(4*n_div)%y
            
                    barra(8*(n_div)-3)%node(1)%x = coord%no(4*n_div-1)%x
                    barra(8*(n_div)-3)%node(1)%y = coord%no(4*n_div-1)%y
                    barra(8*(n_div)-3)%node(2)%x = coord%no(4*n_div)%x
                    barra(8*(n_div)-3)%node(2)%y = coord%no(4*n_div)%y
                    
                    cycle
                end if
                
                
                barra(4*(n_div+k))%conectividades(1) = i-2
                barra(4*(n_div+k))%conectividades(2) = i
                barra(4*(n_div+k)+3)%conectividades(1) = i
                barra(4*(n_div+k)+3)%conectividades(2) = i+1
                
                barra(4*(n_div+k))%node(1)%x = coord%no(i-2)%x
                barra(4*(n_div+k))%node(1)%y = coord%no(i-2)%y
                barra(4*(n_div+k))%node(2)%x = coord%no(i)%x
                barra(4*(n_div+k))%node(2)%y = coord%no(i)%y
            
                barra(4*(n_div+k)+3)%node(1)%x = coord%no(i)%x
                barra(4*(n_div+k)+3)%node(1)%y = coord%no(i)%y
                barra(4*(n_div+k)+3)%node(2)%x = coord%no(i+1)%x
                barra(4*(n_div+k)+3)%node(2)%y = coord%no(i+1)%y
            
            end do
                
        end if
        
        ! cálculo do comprimento de cada barra através das coordenadas das conectividades
        
        do i = 1, n_barras
            barra(i)%comprimento = sqrt((barra(i)%node(2)%y - barra(i)%node(1)%y)**2 + (barra(i)%node(2)%x - barra(i)%node(1)%x)**2)
        end do
        
    end subroutine
    
        !***************************************************************************************************************************************    
        !subroutine secao_barras(n_cant, cant, n_barras, barra)
        !***************************************************************************************************************************************        
        !integer, intent(in) :: n_cant                                ! número de cantoneiras de abas iguais, de acordo com a tabela da AISC de novembro de 2017
        !type(LLS_var), intent(in), allocatable :: cant(:)            ! variável com os dados geométricos de cantoneiras simples de lados iguais
        !integer, intent(in) :: n_barras                              ! número de barras da estrutura treliçada
        !type(barra_trelica), intent(inout), allocatable :: barra(:)  ! vetor com dados geométricos e matriciais das barras da estrutura
        
        !integer :: i = 0
        
        
        ! agrupamento das barras da estrutura reticulada em banzo inferior, banzo superior, montante e diagonal
          
        !integer, allocatable :: i_banzo_superior(:)                  ! armazena os numeros que indicam barras que pertencem ao banzo superior
        !integer, allocatable :: i_banzo_inferior(:)                  ! armazena os numeros que indicam barras que pertencem ao banzo inferior
        !integer, allocatable :: i_diagonal(:)                        ! armazena os numeros que indicam barras que pertencem à diagonal
        !integer, allocatable :: i_montante(:)                        ! armazena os numeros que indicam barras que pertencem ao montante
        !character(20), allocatable :: nome_cant_trelica(:)
        
        ! allocate(nome_cant_trelica(n_barras))
           
        !if(h1>0) then
            
        !    allocate(i_banzo_superior(2*n_div))
        !    allocate(i_montante(2*n_div+1))
        !    allocate(i_diagonal(2*n_div))
        !   allocate(i_banzo_inferior(2*n_div))
            
        !    i_montante(1) = 1
        !    i_banzo_inferior(1)= 2
        !    i_diagonal(1) = 3
        !    i_banzo_superior(1) = 4
            
        !    do i = 1, 2*n_div-1
                
        !        i_montante(i+1) = i_montante(i) + 4
        !        i_banzo_inferior(i+1) = i_banzo_inferior(i) + 4
        !        i_diagonal(i+1) = i_diagonal(i) + 4
        !        i_banzo_superior(i+1) = i_banzo_superior(i) + 4
            
        !    end do
        !        i_montante(i+1) = i_montante(i) + 4
                
        !        nome_cant_trelica(1:n_barras) = ' '
        !        nome_cant_trelica(i_banzo_superior) = cant(20)%name
        !        nome_cant_trelica(i_banzo_inferior) = cant(20)%name
        !        nome_cant_trelica(i_diagonal) = cant(20)%name
        !        nome_cant_trelica(i_montante) = cant(20)%name
                
        !       do i = 1, n_barras
        !            barra(i)%s%secao = LLS_propriedades_cantoneira_lista(n_cant, cant, nome_cant_trelica(i))
        !            call LLS_propriedades_geometricas (barra(i)%s%secao)
         !           call LLSt_propriedades_geometricas (barra(i)%s)
         !       end do
        !else
                
                
        !end if
        !end subroutine
        
        !********************************************************************************************************************
        subroutine carga_pp_barras (rho, h1, n_div, num_nos, n_barras, barra, cond_cont)
        !********************************************************************************************************************
        real(8), intent(in) :: rho                                      ! massa específica do material
        real(8), intent(in) :: h1                                       ! altura do montante de extremidade da cobertura treliçada
        integer, intent(in) :: n_div                                    ! nº de divisões de L/2
        integer, intent(in) :: num_nos                                  ! número total de nós da estrutura treliçada
        integer, intent(in) :: n_barras                                 ! número total de barras da treliça
        type(barra_trelica), intent(inout), allocatable :: barra(:)     ! vetor com dados geométricos e matriciais das barras da estrutura
        type(cond_contorno), intent(out), allocatable :: cond_cont(:)   ! vetor com as condições de contorno de cada nó da estrutura
        
        ! Variaveis internas
        real(8), allocatable :: peso_regiao(:)
        integer :: i=0, n=0
        
        allocate(peso_regiao(4*n_div))
        allocate(cond_cont(num_nos))
        
         
            ! Zerando variáveis
            do i = 1, num_nos
                cond_cont(i)%carga_vento(:) = 0
                cond_cont(i)%carga_telha_telhado(:) = 0
                cond_cont(i)%carga_sobrecarga(:) = 0
                cond_cont(i)%carga_pp(:) = 0
                cond_cont(i)%glc(:) = 0
            end do
            
            ! Atribuição de peso a cada barra
            do i = 1, n_barras
                if(barra(i)%tipo == "2L" .OR. barra(i)%tipo == "2L_cruz") then
                    barra(i)%peso = rho*barra(i)%s%A*barra(i)%comprimento
                else
                    barra(i)%peso = rho*barra(i)%s%secao%A*barra(i)%comprimento
                end if
            end do
            
            do i = 1, n_barras
                cond_cont((barra(i)%conectividades(1)))%carga_pp(2) = cond_cont((barra(i)%conectividades(1)))%carga_pp(2) -barra(i)%peso/2
                cond_cont((barra(i)%conectividades(2)))%carga_pp(2) = cond_cont((barra(i)%conectividades(2)))%carga_pp(2) -barra(i)%peso/2
            end do
            
            ! cálculo do peso da estrutura, dividida em 4*n_div regiões
                !do i = 1, n_div
                 !   peso_regiao(2*i-1) =   rho*SUM(barra(1 +n*4: 1+n*4+3)%comprimento*barra(1 +n*4: 1+n*4+3)%s%A)/2
                 !  peso_regiao(2*i) = rho*SUM(barra((2)+n*4: (2)+n*4+3)%comprimento*barra((2)+n*4: (2)+n*4+3)%s%A)/2
                 !  peso_regiao(4*n_div-2*(i-1)) = peso_regiao(2*i-1)
                 !   peso_regiao(4*n_div-2*(i-1)-1) = peso_regiao(2*i)
                !    n=n+1
                !end do
               ! 
                !do i = 1, n_div+1
               !     if ( i == 1) then
                !        cond_cont(2*i)%carga(2) = -peso_regiao(1)
                !        cond_cont(2*(2*n_div+1))%carga(2) = -peso_regiao(1)
                !!        cycle
                 !   end if
                 !   cond_cont(2*i)%carga(2) = -(peso_regiao(2*(i-1)) + peso_regiao(2*(i-1)+1))
                 !   cond_cont(2*(2*n_div+2) - 2*i)%carga(2) = cond_cont(2*i)%carga(2)
               ! end do
                
            
            ! apoio da esquerda: apoio DUPLO / apoio da direita: apoio SIMPLES
            cond_cont(1)%glc(:) = 1
            cond_cont(num_nos-1)%glc(2) = 1
        
        end subroutine
        
        !********************************************************************************************************************
        subroutine lista_tercas (unidade, terca)
        !********************************************************************************************************************
        
        integer, intent(in) :: unidade                          ! unidade correspondente ao arquivo onde os dados das terças estão armazenados
        type(terca_list), intent(out), allocatable :: terca(:)  ! Variavel que armazena as tercas e suas características identificadoras
        
        ! Variaveis internas
        integer :: stat, num_tercas
        character(20) :: junk
        integer :: i=0, n=0
        
        do
            read(unidade, iostat=stat) junk
            num_tercas=stat-1
            if (stat /= 0) exit
        end do        
        
        allocate(terca(num_tercas))
        rewind(unidade)
        read(unidade,*) junk
        
        do i = 1, num_tercas
            read(unidade, *) terca(i)%nome, terca(i)%carga_linear, terca(i)%A
            terca(i)%carga_linear = terca(i)%carga_linear/10000
        end do
        
        
        end subroutine
        
        !********************************************************************************************************************
        subroutine carga_tercas_telhas_barras (h1, dist_trelica, n_div, num_nos, barra, terca, cond_cont)
        !********************************************************************************************************************
        real(8), intent(in) :: h1                                         ! altura do montante de extremidade da cobertura treliçada
        real(8), intent(in) :: dist_trelica                               ! distância do vão entre dois pórticos consecutivos [cm]
        integer, intent(in) :: n_div                                      ! nº de divisões de L/2
        integer, intent(in) :: num_nos                                    ! número total de nós da estrutura treliçada
        type(barra_trelica), intent(in), allocatable :: barra(:)          ! vetor com dados geométricos e matriciais das barras da estrutura
        type(terca_list), intent(in), allocatable :: terca(:)             ! Variavel que armazena as tercas e suas características identificadoras
        type(cond_contorno), intent(inout), allocatable :: cond_cont(:)   ! vetor com as condições de contorno de cada nó da estrutura
        
        ! Variaveis internas -----------------
        real(8) :: carga_terca                      ! carga das terças (KN/cm)
        real(8) :: carga_telha = 0.0000120d0       ! carga das telhas (KN/cm²)
        real(8) :: d_beiral =  50.d0                ! dimensão do beiral paralela ao banzo superior (cm)
        real(8) :: carga_beiral
        real(8), allocatable :: regiao_telha(:)               
        real(8) :: area_telha
        integer :: i=0, n=0
        

        
        !carga das tercas --------------------
        carga_terca = terca(1)%carga_linear * dist_trelica
        
        do i=1, num_nos/2
            if (i==1 .OR. i==num_nos/2) then
                cond_cont(2*i)%carga_telha_telhado(2) = - carga_terca/2
                cycle
            else if (2*i == 2*(n_div+1)) then
                cond_cont(2*i)%carga_telha_telhado(2) = - carga_terca*2
                cycle
            end if
            cond_cont(2*i)%carga_telha_telhado(2) =  - carga_terca
        end do       
        
        !carga das telhas -----------------------
        allocate(regiao_telha(num_nos/2))
        
        area_telha = L/(2*n_div*cos(theta*pi/180)) * dist_trelica
        carga_telha = carga_telha*area_telha
        carga_beiral = d_beiral*dist_trelica*0.0000120d0
        
        do i = 1, num_nos/2
            if(i == 1 .OR. i == num_nos/2) then
                cond_cont(2*i)%carga_telha_telhado(2) = cond_cont(2*i)%carga_telha_telhado(2) - carga_telha/2 - carga_beiral
                cycle
            end if
            cond_cont(2*i)%carga_telha_telhado(2) = cond_cont(2*i)%carga_telha_telhado(2) - carga_telha
            
        end do
        carga_telha = 0.d0
        
        end subroutine
        
        
        character(len=20) function str(k)
!       "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
        end function str
        
        
        
    end module Estrutura_Trelica