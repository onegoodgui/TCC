 module Analise_trelica_trapezoidal

    use Matriz_Rigidez
    use CargaVento
    use Estrutura_Trelica
    
    implicit none
    
    ! Variaveis ---------------------------------------------------------------------------------------------------------------------
    integer :: iii
    integer, parameter :: uni_dad_estrut = 150, uni_entrada = 300, uni_dad_tercas = 500
    !Variaveis LLS -----------------------------------------------------------------------------------------------------------------
    integer :: n_cant
    type(LLS_var), allocatable :: cant(:)
    !Variaveis barra ---------------------------------------------------------------------------------------------------------------
    integer :: n_barras
    !integer :: num_nos
    real(8) :: rho = 0.000078d0                     ! Massa espec�fica do a�o [KN/cm�]
    real(8) :: fy = 0.25d0                          ! tens�o de escoamento do a�o [KN/mm^2)
    !real(8) :: E = 200.d0                           ! M�dulo de elasticidade do a�o
    real(8) :: G = 77.d0                            ! modulo de elasticidade transversal do a�o [KN/mm^2]
    ! Variables S1-------------------------------------------------------------------------------------------------------------------
    character(4) :: caso                            ! a) Terreno plano ou fracamente acidentado; b) Taludes e morros; c) Vales profundos
    !integer :: theta                               ! Inclina��o m�dia do talude ou enconsta de morro
    real(8) :: z                                    ! Altura medida a partir da superf�cie do terreno no ponto considerado
    real(8) :: d                                    ! diferen�a de n�vel entre a base e o topo do talude ou morro
    namelist /Fator_Top_S1/ caso, theta, z, d       ! namelist com dados de entrada para avaliar o fator  topogr�fico S1
    ! Variables S2-------------------------------------------------------------------------------------------------------------------   
    character :: cat                                ! categoria referente ao fator S2
    character :: cls                                ! classe referente ao fator S2
    namelist /Fator_Top_S2/ cat, cls                ! namelist com dados de entrada para avaliar o fator  topogr�fico S2
    type(tabela_S2) :: tab_S2(5)                    ! vetor de armazenamento dos dados da tabela num�rica do fator S2
     ! Variables press�o din�mica----------------------------------------------------------------------------------------------------
    real(8) :: Vo                                   ! Velocidade caracter�stica do vento da regi�o, em m/s
    namelist /pressao_dinamica/ Vo                  ! namelist com os dados de entrada para avaliar a press�o din�mica
    ! Variables Coef. press�o--------------------------------------------------------------------------------------------------------
    real(8) :: h                                    ! comprimento da altura da coluna de sustenta��o + montante imediatamente acima / vista frontal
    real(8) :: b                                    ! menor dimens�o da estrutura em vista superior
    real(8) :: a                                    ! maior dimens�o da estrutura em vista superior
    integer :: ang_cobertura                        ! inclina��o da viga de cobertura treli�ada entre montante externo e intermedi�rio
    namelist /Cofc_pressao/ h, a, b, ang_cobertura  ! namelist com os valores de entrada para avaliar Ce e Cpe
    !Variables Press�o Vento---------------------------------------------------------------------------------------------------------
    type(Coef_pressao) :: pressao_vento(3)            ! Vari�vel de armazenamento dos coeficientes Ce e Cpe
    real(8) :: coef_min_succao(2), coef_max_succao(2) ! Ce + Ci
    real(8) :: F(2)                                   ! Carga distribuida do vento
    !Variables Dados Geom�tricos ----------------------------------------------------------------------------------------------------
    namelist /Dados_geometricos/ L, h1, dist_trelica, n_div    ! namelist com os valores de entrada referentes �s imposi��es geom�tricas da estrutura
    !Variables Combina��o de a��es --------------------------------------------------------------------------------------------------
    type(comb_acoes) :: comb_acao(5)
    !Variables Otimiza��o
    
    
    integer, parameter :: uni_dad_vento = 100                     ! nome dado para n�mero que identifica arquivo externo
    integer :: i                                                  ! variavel de recurso c�clico
    integer :: n                                                  ! vari�vel de recurso c�clico 
    integer :: junk1
    character :: junk2
    
    contains
    
    !****************************************************************************************************************************************
    subroutine analise_resistencia_axial(barra, fy, E, G, forca_axial, Nrd)
    !****************************************************************************************************************************************
    
    type(barra_trelica), intent(inout), allocatable :: barra(:)     ! Se��o das barras (aplicando simetria), conforme a numera��o e especifica��o da tabela da AISC
    real(8), intent(in) :: fy                                       ! Tens�o de escoamento do a�o [KN/mm^2]
    real(8), intent(in) :: E                                        ! Modulo de elasticidade do a�o [KN/mm^2]
    real(8), intent(in) :: G                                        ! Modulo de elasticidade transversal do a�o [KN/mm^2]
    real(8), intent(in), allocatable :: forca_axial(:)              ! Esfor�o axial nas barras de treli�a [KN]
    real(8) :: Ncrd(n_barras)                                       ! For�a de compress�o resistente de calculo das barras de treli�a [KN]
    real(8) :: Ntrd(n_barras)                                       ! For�a de tra��o resistente de calculo das barras de treli�a [KN]
    real(8), intent(out), allocatable :: Nrd(:,:)                   ! For�a resistente de calculo das barras de treli�a [KN] 

    
    ! Vari�veis internas ------------------
    integer :: i = 0
    integer :: ii = 0
    
    
    allocate(Nrd(6, n_barras))
    
 ! trecho do programa que calcula os esfor�os normais de compress�o e tra��o para a carga externa descendente
 !***********************************************************************************************************
    do i = 1, 6
        do n = 1, n_barras
            if(forca_axial(i)>=0) then
                Ntrd(n) =  barra(n)%s%A*fy/1.1d0
                Nrd(i,n) = Ntrd(n)
            else
                call LLSt_NcRd (barra(n)%s, E, G, fy, barra(n)%comprimento, barra(n)%comprimento, barra(n)%comprimento, NcRd(n))
                Nrd(i,n) = Ncrd(n)
            end if
        end do
    end do
 
    end subroutine
    
    !**********************************************************************************************************************************
    subroutine funcao_objetivo(nd, nc, Vd, Vc, fob)
    !**********************************************************************************************************************************
        integer, intent(in) :: nd                          ! numero de variaveis diescretas
        integer, intent(in) :: nc                          ! numero de variaveis continuas
        integer, intent(in) :: Vd(nd)                      ! valor das variaveis discretas
        real(8), intent(in) :: Vc(nc)                      ! valor das variaveis continuas
        real(8) :: fob                                     ! fun��o objetivo
        real(8), allocatable :: Nrd(:,:)                   ! For�a resistente de calculo das barras de treli�a [KN] 
    
        !Atribui coordenadas para os n�s da estrutura
            call node_barras(n_div, h1, num_nos, coord)
            
        !Estabelece o n�mero de barras e suas conectividades
            call barras_trelica (n_div, h1, coord, n_barras, barra)
            
        !Atribui uma se��o transversal para cada barra de acordo com a lista de se��es da AISC
            call secao_barras(n_cant, cant, n_barras, barra)
    
        !Calcula a a magnitude de carga de peso pr�prio em cada n� da estrutura
            call carga_pp_barras (rho, h1, n_div, num_nos, n_barras, barra, cond_cont)

        !Calcula a carga das ter�as e telhas nos n�s superiores
            call carga_tercas_telhas_barras (h1, dist_trelica, n_div, num_nos, barra, terca, cond_cont)
    
        !Armazena os coeficientes Ce e Cpe da NBR 6123:2013
            call coeficiente_pressao (pressao_vento)
    
        !Combina os valores de Ce e Ci e guarda as combina��es com o valor m�x e min de succ�o
            call combinacao_carregamento_vento (h, a, b, ang_cobertura, pressao_vento, coef_max_succao, coef_min_succao)
    
        !Calcula a carga distribu�da na cobertura devido � a��o do vento
           ! select case(coeficiente_vento)
           !     case(0)
           !         call carga_distribuida_vento(Vo, S1, S2, coef_max_succao, barra, F)
           !     case(1)
           !         call carga_distribuida_vento(Vo, S1, S2, coef_min_succao, barra, F)
           ! end case
        !Calcula as cargas nodais da a��o do vento na estrutura
            call carga_vento_sobrecarga_nos(theta, barra, F, cond_cont)
    
        !Calcula as combina��es de a��es para Combina��es �ltimas Normais 
            call combina_acoes (cond_cont, comb_acao)
    
        !Matriz de rigidez - vetor de deslocamentos fixos e livres / vetor de cargas
            !if (coeficiente_vento == 0) then
            !call vetor_Y_vglc (comb_, Y, v_glc)
            
            allocate(MatRigid(ngln*num_nos,ngln*num_nos))
            MatRigid = 0.d0
    
            do i =1, n_barras
                call matriz_rigidez_local(barra(i)%s%A, E, barra(i)%comprimento, k_local)
                call matriz_trasformacao(barra(i)%node(1)%x,barra(i)%node(1)%y, barra(i)%node(2)%x,barra(i)%node(2)%y, barra(i)%comprimento, T)
                call matriz_rigidez_global(k_local, T, k)
                call matriz_rigidez_estrutura(barra(i)%conectividades(1), barra(i)%conectividades(2), k, MatRigid)
            end do 
    
            call matriz_rigidez_nos_restritos (nglc ,v_glc, MatRigid)
            call eliminacao_gauss(num_nos, MatRigid, Y)
            call deslocamentos_matriz_rigidez(nglc, v_glc, MatRigid, Y, deslocamentos)
            call forca_axial_barras(barra, n_barras, deslocamentos, forca_axial)
            call analise_resistencia_axial(barra, fy, E, G, forca_axial, Nrd)
            
    end subroutine
    
    end module