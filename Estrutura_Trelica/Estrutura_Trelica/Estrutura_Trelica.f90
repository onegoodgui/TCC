
    module Estrutura_Trelica

    ! módulo para calculo dos esforços numa estrutura de treliça trapezoidal/triangular, de número de
    implicit none

    
    ! Variáveis TYPE -------------------------------------------------------------------------------------
    type coordenadas
        real(8) :: x
        real(8) :: y
    end type
    
    type node
        type(coordenadas), allocatable :: no(:)
    end type
    
    type barra_trelica
        type(coordenadas) :: node(2)
        real(8) :: conectividades(2)
        real(8) :: comprimento
    end type
        
    
    type estrutura
        real(8), allocatable :: montante(:)
        real(8), allocatable :: banzo_inferior(:)
        real(8), allocatable :: diagonal(:)
        real(8), allocatable :: banzo_superior(:)
    end type
    
    ! Variaveis para montagem da estrutura ---------------------------------------------------------------
    real(8) :: L = 18                                        ! vão horizontal entre as extremidades do pórtico
    real(8) :: h1 = 0                                       ! altura do montante de extremidade
    integer :: n_div = 3                                    ! divisões do comprimento referente a L/2
    integer :: theta =15                                    ! inclinação do banzo superior
    !namelist /dad_entrada_trelica/ L, h1, n_div, theta   ! namelist com os dados necessários para montagem da estrutura
    type(node) :: coord
    type(barra_trelica), allocatable :: barra(:)
    
  
    
    
    real(8), parameter :: pi = 4.d0*atan(1.d0)
    
    contains
    
    !subroutine geometria_estrutura (l, h1, theta, n_div, estrutura)
    !
    !real(8), intent(in) :: l                                       ! vão horizontal entre as extremidades do pórtico
    !real(8), intent(in) :: h1                                      ! altura do montante de extremidade
    !real(8), intent(in) :: theta                                   ! inclinação do banzo superior
    !real(8), intent(in) :: n_div                                   ! divisões do comprimento referente a l/2
    !type(estrutura), intent(out) :: trelica                        ! estrutura da trelica, dividida em banzo inferior, banzo superior, montante e diagonal
    !
    !allocate(trelica%banzo_inferior(n_div))
    !allocate(trelica%banzo_superior(n_div))
    !allocate(trelica%diagonal(n_div))
    !
    !allocate(trelica%banzo_inferior(n_div))
    
    
    !Subrotina que numera os nós da treliça e calcula as coordenadas de cada nó na sua respectiva numeração
    subroutine node_barras(n_div, h1, coord)
    
    integer(4), intent(in) :: n_div
    real(8), intent(in) :: h1
    type(node), intent(out) :: coord
    
    integer :: i = 0, n = 0
    
    
        if (h1>0) then
            
            !---------COBERTURA TRAPEZOIDAL -----------
            allocate(coord%no(4*n_div + 2))
             
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
                coord%no(2*(i))%y = h1 + (i-1)*tan(theta*pi/180)
             end do
        else
            !---------COBERTURA TRIANGULAR -----------
            allocate(coord%no(4*n_div))
            
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
    
    
    subroutine barras_trelica (n_div, h1, coord, barra)
    
        integer(4), intent(in) :: n_div                                ! nº de divisões de L/2
        real(8), intent(in) :: h1                                      ! altura do montante de extremidade da cobertura
        type(node), intent(in) :: coord                                ! vetor com os nós e suas respectivas coordenadas
        type(barra_trelica), intent(out), allocatable :: barra(:)      ! vetor com as barras, suas conectividades e seu comprimento
        
        integer :: n=0, i=0, k=0                                           ! contadores
        
        if(h1>0) then
            !---------COBERTURA TRAPEZOIDAL -----------
            allocate(barra(8*n_div+1))
            
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
            
            barra(4*n_div+1)%conectividades(1) = i+1
            barra(4*n_div+1)%conectividades(2) = i+2
            
            barra(4*n_div+1)%node(1)%x = coord%no(i+1)%x
            barra(4*n_div+1)%node(1)%y = coord%no(i+1)%y
            
            barra(4*n_div+1)%node(2)%x = coord%no(i+2)%x
            barra(4*n_div+1)%node(2)%y = coord%no(i+2)%y
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
            
        else
            !---------COBERTURA TRIANGULAR -----------
            allocate(barra(8*n_div-3))
            
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
        
        
        
        
        
        
    end subroutine
    
        
    end module Estrutura_Trelica

