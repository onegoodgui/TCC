module LLSt
    
    use LLS
    
implicit none
    
    !*****************************************************************
    !m�dulo para dimensionamento de elementos compostos por cantoneiras de a�o de abas iguais, com se��o transversal em T. 
    !A se��o de refer�ncia possui como eixo x o eixo colinear as abas horizontais alinhadas, e como eixo y, o eixo situado a meia distancia das abas verticais
    !*****************************************************************
    
    
    !*****************************************************************
    !declara��o da vari�vel tipo
    !*****************************************************************
    type LLSt_var
        
    !*****************************************************************
    !    Entradas
    !*****************************************************************
    type(LLS_var) :: secao        ! propriedades geometricas da secao cantoneira simples 
    real(8) :: h = 0.d0           ! Altura da aba da cantoneira  
    real(8) :: dist =0.d0         ! distancia entre as faces verticais das se�oes cantoneiras
    real(8) :: t_o = 0.d0         ! Espessura da aba
    real(8) :: A = 0.d0           ! �rea da se��o T
    real(8) :: x_g = 0.d0         ! para os eixos de referencia, x_g =0
    real(8) :: y_g = 0.d0         ! dist�ncia entre a face da aba horizontal e o eixo x da cantoneira simples
    real(8) :: r_max = 0.d0       ! Raio de giracao maximo
    real(8) :: r_min = 0.d0       ! Raio de giracao minimo
    real(8) :: Ix = 0.d0          ! Momento de in�rcia em rela��o ao eixo x
    real(8) :: Iy = 0.d0          ! Momento de in�rcia em rela��o ao eixo y
    real(8) :: J = 0.d0           ! Momento de in�rcia a tor��o pura
    real(8) :: Cw = 0.d0          ! Constante de empenamento
    real(8) :: r_o= 0.d0          ! raio de giracao polar em rela��o ao centro de cisalhamento
    
    !******************************************************************
    !  Componentes geom�tricos da se��o
    !******************************************************************
    
    
    real(8) :: Ixy = 0.d0         ! Produto de in�rcia em rela��o aos eixos x e y
    real(8) :: rx = 0.d0          ! Raio de gira��o em rela��o ao eixo x
    real(8) :: ry = 0.d0          ! Raio de gira��o em rela��o em eixo y
    real(8) :: Wx_s = 0.d0        ! M�dulo resistente a flex�o em rela��o ao eixo x superior
    real(8) :: Wx_i = 0.d0        ! M�dulo resistente a flex�o em rela��o ao eixo x inferior
    real(8) :: Imax = 0.d0        ! Momento de in�rcia central m�ximo
    real(8) :: y_o = 0.d0         !distancia entre o centroide e o centro de cisalhamento da se��o transversal da cantoneira (y_o)
    real(8) :: lambda = 0.d0       !  indice de esbeltez da aba da cantoneira
    
    
    
    end type LLSt_var
    !************************************************************************

    
    
    !************************************************************************
    !Declara��o de par�metros
    !************************************************************************
    
    real(8), parameter, private:: pi = 4.d0*atan(1.d0)    !N�mero, se uma circunfer�ncia tem per�metro p e diametro d, pi = p/d
    real(8), parameter, private :: gama_a1 = 1.10         ! coeficiente de redu�ao de resistencia para combinacao normal de a�oes
    
    !************************************************************************
    
    
    
    
    contains
    

    
    !******************************************************************************************************
    !******************************************************************************************************
    !******************************************************************************************************
    !******************************************************************************************************
       
    !********************************************************************************************************************************
    ! subrotina que armazena propriedades geometricas derivadas da secao simples e calcula propriedades exclusivas da configuracao "T"
    !********************************************************************************************************************************
    subroutine LLSt_propriedades_geometricas (cant_t)
   

    type(LLSt_var),intent(inout) :: cant_t
    
 
       
 cant_t%A = 2*cant_t%secao%A
 cant_t%h = cant_t%secao%h
 cant_t%t_o = cant_t%secao%t_o
 cant_t%x_g = 0
 cant_t%y_g = cant_t%secao%y_g
 cant_t%Ix = 2*cant_t%secao%Ix         
 cant_t%rx = cant_t%secao%rx           
 cant_t%r_min = cant_t%rx                                           
 cant_t%J = 2*cant_t%secao%J
 cant_t%Cw = 2*cant_t%secao%Cw

 


       
    !! subrotina para calculo da momento de inercia em rela��o ao eixo y da se��o composta (Iy)
    cant_t%Iy = 2.d0*cant_t%secao%Iy + 2.d0*cant_t%secao%A*(cant_t%secao%x_g + cant_t%dist/2.d0)**2
    !
    !
    !! c�lculo do momento de inercia maximo (Imax)
     cant_t%Imax = cant_t%Iy
    !
    !!c�lculo do raio de giracao maximo (r_max)
    cant_t%r_max = sqrt(cant_t%Imax/cant_t%A)
    cant_t%ry = cant_t%r_max
    !
    ! ! c�lculo do m�dulo resistente a flex�o superior em rela��o ao eixo x (Wx_s)
      cant_t%Wx_s = cant_t%Ix/(cant_t%h-cant_t%y_g)
    
    !
    ! ! c�lculo do m�dulo resistente a flex�o inferior em rela��o ao eixo x (Wx_i)
      cant_t%Wx_i = cant_t%Ix/cant_t%y_g
    
    !
    ! ! c�lculo do momento de in�rcia a tor��o pura (J)
    !cant_t%J = (4.d0/3.d0)*(cant_t%h - cant_t%t_o/2.d0)*cant_t%t_o**3
    
    ! 
    ! ! c�lculo da constante de empenamento (Cw)
    ! cant_t%Cw = (1.d0/9.d0)*((cant_t%h - cant_t%t_o/2.d0)**3.d0)*(cant_t%t_o**3.d0)
     
    !
    ! !c�lculo da distancia entre o centroide e o centro de cisalhamento da se��o transversal da cantoneira (y_o)
      cant_t%y_o = cant_t%secao%y_g - cant_t%secao%t_o/2
     
    ! 
    ! !c�lculo do raio de giracao polar em rela��o ao centro de cisalhamento (r_o)
      cant_t%r_o = sqrt(cant_t%ry**2.d0 + cant_t%rx**2.d0 + cant_t%y_o**2 )
    ! 
    ! ! c�lculo da esbeltez da aba (lamba) 
      cant_t%lambda = cant_t%h/cant_t%t_o
     

   
end subroutine
    
     !******************************************************************************************************
     !******************************************************************************************************
    
    
    
    
    
    
    
    
    
     !******************************************************************************************************
     ! c�lculo do fator de reducao para flambagem local
     !******************************************************************************************************
     
     
     subroutine LLSt_Qs(cant_t, E, fy, Qs) ! c�lculo do fator de reducao para flambagem local
     
     type(LLSt_var), intent(in) :: cant_t           
     real (8), intent(in) :: E                       ! M�dulo de elasticidade do a�o
     real (8), intent(in) :: fy                      ! tens�o de escoamento do a�o (conforno tabela A1.4 da NBR8800)
     real (8), intent(out) :: Qs                     ! coeficiente de reducao aplicavel a se�oes em que os componentes tem indice de esbeltez superior aos valores da tabela F.1 da NBR8800
     real(8) :: lambda_r                             ! esbeltez limite (tabela F1 da NBR8800)
     real(8) :: lambda_e                             ! esbeltez limite elastica (tabela F1 da NBR8800)
    
    
    lambda_r = 0.45d0*sqrt(E/fy)                     ! esbeltez limite 'r' da aba da cantoneira. abaixo do qual n�o ocorre flambagem local. Conforme livro Pfeil - Estruturas de Aco, Dimensionamento Pr�tico; pag 133
    lambda_e = 0.91d0*sqrt(E/fy)                     ! esbeltez limite 'e' da aba da cantoneira. acima do qual ocorre flambagem local inel�stica. Conforme livro Pfeil - Estruturas de Aco, Dimensionamento Pr�tico; pag 133

    
     if( cant_t%lambda <= lambda_r ) then 
         
         Qs = 1.0d0 ! flambagem local inexistente. apenas flambagem global
         
       
         
     end if
     
     if (cant_t%lambda <= lambda_e .AND. cant_t%lambda > lambda_r) then
     
     Qs = 1.340d0 - 0.76d0*(cant_t%h/cant_t%t_o)*sqrt(fy/E)  
 
     
     else if (cant_t%lambda > lambda_e) then
         
     Qs = (0.53d0*E)/(fy*((cant_t%h/cant_t%t_o)**2))
 
     
     end if
     


     end subroutine
!     
!     !******************************************************************************************************
!     !******************************************************************************************************
!     
!     
!     
!     
     !******************************************************************************************************
     ! calculo da for�a de flambagem cr�tica
     !******************************************************************************************************
     
     
    subroutine LLSt_Ncr(cant_t, E, G, KxLx, KyLy, KzLz, Ncr) ! calculo da for�a de flambagem cr�tica
    
    type(LLSt_var),intent(in) :: cant_t
                                                        
    real(8), intent(in) :: G                                ! modulo de elasticidade transversal do a�o (KN/cm2)
    real(8),intent(in):: E                              ! m�dulo de elasticidade do a�o 
    real(8),intent(in):: KxLx                           ! comprimento de flambagem em rela��o ao eixo principal central x
    real(8),intent(in):: KyLy                           ! comprimento de flambagem em rela��o ao eixo principal central y 
    real(8),intent(in):: KzLz                           ! comprimento de flambagem em rela��o ao eixo z
    real(8) :: Nex                                      ! for�a axial de flambagem elastica por flexao em rela��o ao eixo central principal x
    real(8) ::   Ney                                    ! for�a axial de flambagem elastica por flexao em rela��o ao eixo central principal y
    real(8) ::   Nez                                    ! for�a axial de flambagem elastica por tor��o em rela��o ao eixo z
    real(8) ::  Neyz                                    ! for��o de flambagem el�stica por flexo-tor��o
    real(8), intent(out) :: Ncr                         ! for�a de flambagem critica
 
    
    
    
     Nex = ((pi**2)*E*cant_t%A*cant_t%r_min**2)/(KxLx**2)
     
     Ney = ((pi**2)*E*cant_t%A*cant_t%r_max**2)/(KyLy**2)
    
     Nez = (1.d0/cant_t%r_o**2)*(G*cant_t%J + (pi**2)*(E*cant_t%Cw)/(KzLz**2))
     
     Neyz = ((Ney + Nez)/(2.d0*(1.d0 - (cant_t%y_o/cant_t%r_o)**2)))*(1.d0 - sqrt(1.d0 - (4.d0*Ney*Nez*(1.d0 - (cant_t%y_o/cant_t%r_o)**2))/(Ney + Nez)**2))
     
     Ncr = min(Nex,Neyz)
     

    
    end subroutine
    
    !******************************************************************************************************
    !******************************************************************************************************
    
    
    
    !******************************************************************************************************
    ! calculo da for�a axial de compress�o resistente de c�lculo
    !******************************************************************************************************
    
    subroutine LLSt_NcRd (cant_t, E, G, fy, KxLx, KyLy, KzLz, NcRd)
    
    type(LLSt_var), intent(inout) :: cant_t
    
    real(8), intent(in) :: E                                              ! M�dulo de elasticidade do a�o
    real(8), intent(in) :: G                                              ! modulo de elasticidade transversal do a�o    
    real(8), intent(in) :: fy                                             ! tens�o de escoamento do a�o (conforno tabela A1.4 da NBR8800)
    real(8):: KxLx                                                        ! comprimento de flambagem em rela��o ao eixo principal central x
    real(8):: KyLy                                                        ! comprimento de flambagem em rela��o ao eixo principal central y
    real(8):: KzLz                                                        ! comprimento de flambagem em rela��o ao eixo z
    real(8) :: lambda_0                                                   ! indice de esbeltez reduzido
    real(8) :: qui                                                          ! parametro adimensional que relaciona tens�o cr�tica e tens�o de escoamento do a�o
    real(8), intent(out) :: NcRd                                          ! for�a axial de compress�o resistente de c�lculo
    
    real(8) ::  Qs
    real(8) :: Ncr
    
    !call  LLSt_propriedades_geometricas (cant_t)
    call  LLSt_Qs(cant_t, E, fy, Qs)
    call  LLSt_Ncr(cant_t, E, G, KxLx, KyLy, KzLz, Ncr)
    

    lambda_0 = sqrt((Qs*cant_t%A*fy)/Ncr)
    
    if(lambda_0 <= 1.50d0) then
        
        qui = 0.658d0**(lambda_0**2)
        
    else
        qui = 0.877d0/(lambda_0)**2
         
    end if
    
       NcRd = (Qs*cant_t%A*qui*fy)/gama_a1
     
    
    end subroutine
    

    
    
    
    
        end module LLSt