module Matriz_Rigidez
    
use LLSt
use LLS
use Estrutura_Trelica    

implicit none


real(8) :: E = 20000.d0                         ! Módulo de elasticidade do aço [KN/cm^2]
!real(8) :: L                                    ! comprimento da barra analisada
integer :: i_inicio                             ! indice do ponto inicial
integer :: j                                    ! indice do ponto final
real(8) :: xi                                   ! coordenada x do nó inicial da barra
real(8) :: yi                                   ! coordenada y do nó inicial da barra
real(8) :: xf                                   ! coordenada x do nó final da barra
real(8) :: yf                                   ! coordenada y do nó final da barra
real(8) :: cj(2)                                ! vetor das coordenadas globais em relação ao eixo y do ponto inicial
real(8) :: k_local(4,4)                         ! Matriz de rigidez local
real(8) :: T (4, 4)                             ! Matriz de transformação
real(8) :: k(4, 4)                              ! Matriz de rigidez global
real(8) :: contador
!real(8) :: Vd(4)                                ! variaveis discretas (referentes as dimensoes das seçoes consideradas nos calculos)

!Variaveis Matriz de Rigidez
    integer :: ngln = 2                             ! numero de graus de liberdade por nós
    !integer :: num_nos                             ! numero de nós da estrutura
    integer :: MDC                                  ! maxima diferença na conectividade da malha
    integer :: nne = 2                              ! número de nos por elemento
    integer :: nglc = 3                             ! numero de graus de liberdade com condição de contorno
    integer :: ntgl                                 ! numero total de graus de liberdade
    integer :: nncc                                 ! numero de nos com condição de contorno
    integer :: nnf                                  ! numero de nos com forças atuantes
    integer :: err
    integer, allocatable :: v_glc(:)
    !type(MDSB_pbMEF) :: MDSB_trelica
    real(8), allocatable :: Y(:,:)
    real(8), allocatable :: MatRigid(:,:)
    !namelist/nos_graus_liberdade/ nglc, nncc
    !namelist/nos_forcas/ ntgl, nnf

    contains
    
    
    
! **********************************************************************    
! subroutina que monta a matriz de rigidez local do elemento selecionado
! **********************************************************************
    subroutine matriz_rigidez_local(A, E, L, k_local)
    
    real(8), intent(in) :: A                                   ! variavel que contem propriedades geometricas do elemento
    real(8),intent(in) :: E                                   ! Módulo de elasticidade do aço
    real(8), intent(in) :: L                       ! comprimento da barra analisada
    
    real(8), intent(out) :: k_local(4,4)           ! Matriz de rigidez local
    
    
    k_local(1,1) = (A*E)/L
    k_local(1,2) = 0
    k_local(1,3) = -(A*E)/L
    k_local(1,4) = 0
    k_local(2,1) = 0
    k_local(2,2) = 0
    k_local(2,3) = 0
    k_local(2,4) = 0
    k_local(3,1) = -(A*E)/L
    k_local(3,2) = 0
    k_local(3,3) = (A*E)/L
    k_local(3,4) = 0
    k_local(4,1) = 0
    k_local(4,2) = 0
    k_local(4,3) = 0
    k_local(4,4) = 0
    

    
    end subroutine
    
 ! *************************************************************************
 ! subroutina que constroi a matriz de transformação do elemento selecionado  
 ! *************************************************************************
    subroutine matriz_trasformacao(xi, yi, xf, yf, L, T)
    
    real(8) :: xi 
    real(8) :: yi
    real(8) :: xf
    real(8) :: yf
    real(8) :: c                     ! Variavel que armazena o cosseno referente ao angulo entre a direção do elemento e o eixo x
    real(8) :: s                     ! Variavel que armazena o seno referente ao angulo entre a direção do elemento e o eixo x
    
    real(8), intent(out) :: L        ! Comprimento do elemento
    real(8), intent(out) :: T (4, 4) ! Matriz de transformação
    
    L = sqrt((xf-xi)**2 + (yf-yi)**2)
    c = (xf - xi)/L
    s = (yf - yi)/L
    
    T(1, 1) = c
    T(1, 2) = s
    T(1, 3) = 0
    T(1, 4) = 0
    T(2, 1) = -s
    T(2, 2) = c
    T(2, 3) = 0
    T(2, 4) = 0
    T(3, 1) = 0
    T(3, 2) = 0
    T(3, 3) = c
    T(3, 4) = s
    T(4, 1) = 0
    T(4, 2) = 0
    T(4, 3) = -s
    T(4, 4) = c
    
    end subroutine
    
    
 ! *************************************************************************
 ! subroutina que constroi a matriz de rigidez global do elemento selecionado  
 ! *************************************************************************  
    
    subroutine matriz_rigidez_global(k_local, T, k)
    
    real(8), intent(in) :: k_local(4, 4)
    real(8), intent(in) :: T(4, 4)
    
    real(8), intent(out) :: k(4, 4)
    
    k = matmul(matmul(transpose(T),k_local),T)
    
    
    end subroutine
    
    !********************************************************************************************************************
    subroutine matriz_rigidez_estrutura(no_inicial, no_final, k, MatRigid)
    !********************************************************************************************************************
    integer :: no_inicial
    integer :: no_final
    real(8), intent(in) :: k(4,4)
    real(8), intent(inout) :: MatRigid(:,:)
    
    !Variaveis internas
    integer :: i=0, n=0, j=0
    integer :: conect(4)
    
    conect(1) = 2*no_inicial-1
    conect(2) = 2*no_inicial
    conect(3) = 2*no_final -1
    conect(4) = 2*no_final
    
    do i = 1,4
        MatRigid(conect(:),conect(i)) = MatRigid(conect(:),conect(i)) + k(:,i)
    end do

    end subroutine
    !*************************************************************************************************************
    subroutine vetor_vglc (vetor, v_glc)
    !*************************************************************************************************************
    type(cond_contorno), allocatable, intent(in) :: vetor(:)
    integer, allocatable, intent(inout) :: v_glc(:)
    
    !Variaveis internas
    real(8), allocatable :: aux(:,:)
    integer :: i=0, n=0, k=0, j=0
    
    allocate(aux(num_nos*ngln,1))

    
    do i = 1, num_nos
        do n = 1, ngln
            aux(ngln*i-ngln+n,1) = vetor(i)%glc(n)
            
                if (vetor(i)%glc(n) == 1) then
                k=k+1
                v_glc(k) = ngln*(i-1) + n
                end if
        end do
    end do
    k=0
    end subroutine
    
    !********************************************************************************************************************
    subroutine vetor_Y (num_no,  vetor, carga_pp, carga_telha_telhado, carga_sobrecarga, carga_vento, Y)  ! subrotina que armazena os ntgcl de acordo com as restrições de cada um
    !********************************************************************************************************************
    integer, intent(in) :: num_no
    type(cond_contorno), intent(in) :: vetor(num_no)
    real(8), intent(in)  :: carga_pp(2)
    real(8), intent(in)  :: carga_telha_telhado(2)
    real(8), intent(in)  :: carga_sobrecarga(2)
    real(8), intent(in)  :: carga_vento(2)
    real(8), allocatable, intent(inout) :: Y(:,:)
    
    !Variaveis internas
    real(8), allocatable :: aux(:,:)
    integer :: i=0, n=0, k=0, j=0
    
    allocate(aux(num_nos*ngln,1))

    

        do n = 1, ngln
            !aux(ngln*i-ngln+n,1) = vetor(num_no)%glc(n)
            aux(ngln*num_no-ngln+n,1) = carga_pp(n) + carga_telha_telhado(n) + carga_sobrecarga(n) + carga_vento(n)
            Y(ngln*num_no-ngln+n,1) = aux(ngln*num_no-ngln+n,1)
        end do

    end subroutine
    
    !*******************************************************************************************************************************************************************    
    subroutine matriz_rigidez_nos_restritos (nglc ,v_glc, MatRigid) ! subrotina que altera matriz de coeficientes de acordo com as condições de restrições dos nós especificados
    !*******************************************************************************************************************************************************************
    
    integer, intent(in) :: nglc
    integer, allocatable, intent(in) :: v_glc(:)
    real(8), allocatable, intent(inout) :: MatRigid(:,:)
    
    !Variaveis internas
    integer :: i=0, j=0
    
    do i=1, nglc
        !MatRigid(v_glc(i),:) = 0.d0
        MatRigid(:,v_glc(i)) = 0.d0
        !MatRigid(v_glc(i),v_glc(i)) = 1.d0
    end do
    
    
    end subroutine
    
 
    !********************************************************************************************************************
    subroutine forca_axial_barras(barra, n_barras, deslocamentos, forca_axial)
    !********************************************************************************************************************
    type(barra_trelica), intent(in),allocatable :: barra(:)     ! vetor com as barras, suas conectividades e seu comprimento
    integer, intent(in) :: n_barras
    real(8), intent(in),allocatable :: deslocamentos(:,:)       ! matriz que contém valor dos deslocamentos na direção x e y dos nós da treliça
    real(8), intent(OUT), allocatable :: forca_axial(:)         ! força axial em cada barra
    
    ! Variáveis internas
    integer :: i=0 , n=0
    real(8) :: L_def
    real(8) :: delta_L
    allocate(forca_axial(n_barras))
    
    do i = 1, n_barras
    
        L_def = sqrt(((barra(i)%node(2)%x + deslocamentos(2*barra(i)%conectividades(2)-1,1)) - (barra(i)%node(1)%x + deslocamentos(2*barra(i)%conectividades(1)-1,1)))**2 + ((barra(i)%node(2)%y + deslocamentos(2*barra(i)%conectividades(2),1)) - (barra(i)%node(1)%y + deslocamentos(2*barra(i)%conectividades(1),1)))**2)
        delta_L = L_def - barra(i)%comprimento
        if(barra(i)%tipo == "2L" .OR. barra(i)%tipo == "2L_cruz") then
            forca_axial(i) = barra(i)%s%A*delta_L*E/barra(i)%comprimento
        else
            forca_axial(i) = barra(i)%s%secao%A*delta_L*E/barra(i)%comprimento
        end if
        
    end do
    
    end subroutine

!***********************************************************************************************************    
subroutine eliminacao_gauss(num_nos, matriz, Y)
!***********************************************************************************************************
integer, intent(in) :: num_nos
real(8), intent(inout), allocatable :: matriz(:,:)
real(8), intent(inout), allocatable :: Y(:,:)

!Variáveis internas
integer :: i = 0, j=0, h=0, k=0, n=0
integer  :: ntgl

ntgl = num_nos*ngln
!allocate(deslocamentos(ntgl))

do h = 1,ntgl-1
    if(matriz(h,h) == 0.d0) then
        cycle
    end if
    Y(h,1) = Y(h,1)/matriz(h,h)
    matriz(h,:) = matriz(h,:)/matriz(h,h)
    do j = h + 1, ntgl
        Y(j,1) = Y(j,1) - Y(h,1)*matriz(j,h)
        matriz(j, :) = matriz(j, :) - matriz(j,h)*matriz(h, :)
    end do
end do

end subroutine

!*********************************************************************************************************
subroutine deslocamentos_matriz_rigidez(nglc, v_glc, matriz, Y, deslocamentos)
!*********************************************************************************************************
integer, intent(in) :: nglc
integer, intent(in), allocatable :: v_glc(:)
real(8), intent(inout), allocatable :: matriz(:,:)
real(8), intent(inout), allocatable :: Y(:,:)
real(8), intent(inout), allocatable :: deslocamentos(:,:)

!Variaveis internas
integer :: i=0, j=0, k=0, n=0
integer :: ntgl

ntgl = ngln*num_nos
allocate(deslocamentos(ntgl,1))
deslocamentos(:,1) = 1.d0

do i=1, nglc
    deslocamentos(v_glc(i),1) = 0.d0
end do
    
    deslocamentos(ntgl,1) = Y(ntgl,1)/(matriz(ntgl,ntgl))
    
    
do i =2, ntgl
    if(deslocamentos(ntgl-i+1,1)  == 0.d0) then
       cycle
    end if
    n=i-1
    deslocamentos(ntgl-i+1,1) = (Y(ntgl-i+1,1) - SUM(matriz(ntgl-i+1,ntgl-n+1:ntgl)*deslocamentos(ntgl-n+1:ntgl,1))/matriz(ntgl-i+1,ntgl-i+1))
end do    

end subroutine   
      
end module
    