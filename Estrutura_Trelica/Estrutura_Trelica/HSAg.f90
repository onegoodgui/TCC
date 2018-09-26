MODULE HSAg
!************************************************************************************
!* PROP�SITO:																		*
!* Realiza otimiza��o pelo Algoritmo de Busca pela Harmonia (HSA)					*
!* (HSA = Harmony Search algorithm)           										*
!* 																					*
!* Solu��o de problemas de otimiza��o (lineares ou n�o-lineares)					*
!* Vari�veis inteiras ou vari�veis discretas conectadas a inteiras.					*
!* (linked-discrete variables nonlinear programing LD-NLP)							*
!* Vari�veis reais																	*
!* 																					*
!* 																					*
!*	O problema de otimiza��o deve ser do tipo:										*
!*																					*
!*		Minimizar: F(X)																*
!*		tal que X � um vetor com n compnentes x_i									*
!* 		e x_i_L <= x_i <= x_i_U														*
!* 																					*
!*      F(X) � a fun��o objetivo													*
!* 		No caso de otimiza��es com restri��es F(X) deve contemplar fun��es de 		*
!* 		penalidade																	*
!* 																					*
!* Refer�ncias:  																	*
!*		1) Z.W.Geem,J.-H.Kim,G.V.Loganathan,A new heuristic optimization            *
!*		algorithm: harmony search.Simulation 76(2)(2001)60�68.	                    *
!* 		2) A new meta-heuristic algorithm for continuous engineering optimization   *
!* 		 - harmony search theory and practice.                                      *
!* 		 Comput.Methods Appl. Mech. Engrg.194(2005)3902�3933						*
!* 		3) Optimum design of steel frames using a harmony search algorithm.         * 
!*      Struct. Multidisc Optim (2008) 36:393-401                                   *
!* 																					*
!* 																					*
!* 																					*
!*																					*
!* Autor: Felipe Schaedler de Almeida (Almeida, F. S.)	 			
!************************************************************************************
USE numRAND     !M�dulo numRAND com subrotinas auxiliares para gera��o de valores
                !rand�micos l�gicos, inteiros ou reais em intervalos definidos

IMPLICIT NONE   !Todas as vari�veis devem ser declaradas
PRIVATE         !Todas as rotainas s�o privadas a menos da lista abaixo

!-------------Lista de subrotinas e vari�veis p�blicas-------------
public otimiza_HSA
public iniciar_HSA_sai
public terminar_HSA_sai

public imprimir_HSA_sai_hAF_plt


public HSA_pto
public HSA_par
public HSA_sai


!public fnc


!----------Par�metros internos do algoritmo  -----------------
!(usadas para aloca��o est�tica de algumas vari�veis)
integer,parameter:: maxNvD = 50 !<N�m. m�ximo de vari�veis discretas
integer,parameter:: maxNvC = 50 !<N�m. m�ximo de vari�veis cont�nuas
integer,parameter:: ncoef = 5   !<N�m. de coeficientes das fun��es de varia��o do par�metros
integer,parameter:: nctf = 3    !<N�m. de caracteres usados para definir as fun��es 


!-----------Declara��o das vari�veis tipo derivadas--------------

type HSA_pto !<ponto avaliado pelo algoritmo no espa�o de respostas do problema
    integer,allocatable:: vd(:) !<vari�vel discreta
    real(8),allocatable:: vc(:) !<vari�veis cont�nuas
    real(8):: fob = 0.d0        !<Valor da fun��o objetivo
end type HSA_pto

type HSA_par !<Par�metros do HSA
    
    integer:: nvD  = 0          !<n�mero de vari�veis discretas
    integer:: nvC  = 0          !<n�mero de vari�veis cont�nuas
    integer:: Li_vd(maxNvD) = 0 !<Limite inferior das vari�veis discretas
    integer:: Ls_vd(maxNvD) = 0 !<Limite superior das vari�veis discretas
    real(8):: Li_vc(maxNvC) = 0.d0    !<Limite inferior das vari�veis cont�nuas
    real(8):: Ls_vc(maxNvC) = 0.d0    !<Limite superior das vari�veis cont�nuas  
    
    integer:: maxIT = huge(1)   !N�m. m�ximo de itera��es do algoritmo
    integer:: maxAF = huge(1)   !N�m. m�ximo de avalia��es da fun��o objetivo (sem considerar as repeti��es)
    
    logical:: verificar_parametros = .false. !Indica se deve ser realizada a verifica��o dos par�metros previamente � otimiza��o

    character(2):: contador = 'it'      !<define se usa o n�m. de itera��es ('it' e default) ou o n�m de an�lises da fun��o objetivo (af) como contador
    integer:: maxHMS = 0                !<Tamanho m�ximo da HM
    character(nctf):: HMS_tf  = 'cst'   !<Fun��o de varia��o de HM
    real(8):: HMS_c(ncoef) = 0          !<Coeficientes da fun��o de varia��o de HM
    
    !Par�metros para as vari�veis discretas  ---vd---    
    character(nctf):: HMCRvd_tf(maxNvD) = 'cst'   !<Fun��o de varia��o de HMCR para vd
    real(8):: HMCRvd_c(ncoef,maxNvD) = 0.d0       !<Coef. da fun��o de varia��o de HMCR para vd

    character(nctf):: PARvd_tf(maxNvD) = 'cst'    !<Fun��o de varia��o de PAR para vd
    real(8):: PARvd_c(ncoef,maxNvD) = 0.d0        !<Coef. da fun��o de varia��o de PAR para vd
    
    character(nctf):: PAbw_vd_tf(maxNvD) = 'cst'  !<Fun��o de varia��o de bw para vd
    real(8):: PAbw_vd_c(ncoef,maxNvD) = 0         !<Coef. da fun��o de varia��o de bw para vd
    
    !Par�metros para as vari�veis const�nuas  ---vc---    
     
    character(nctf):: HMCRvc_tf(maxNvC) = 'cst'   !<Fun��o de varia��o de HMCR para vc
    real(8):: HMCRvc_c(ncoef,maxNvC) = 0.d0       !<Coef. da fun��o de varia��o de HMCR para vc
    
    character(nctf):: PARvc_tf(maxNvC) = 'cst'    !<Fun��o de varia��o de PAR para vc
    real(8):: PARvc_c(ncoef,maxNvC) = 0.d0        !<Coef. da fun��o de varia��o de PAR para vc
    
    character(nctf):: PAbw_vc_tf(maxNvC) = 'cst'  !<Fun��o de varia��o de bw para vc
    real(8):: PAbw_vc_c(ncoef,maxNvC) = 0.d0      !<Coef. da fun��o de varia��o de bw para vc
    
    !Par�metros para controle de impress�o de dados durante a otimiza��o
    integer:: uni_imp_par = 0   !<Unid. para impress�o de par�metros do HSA (0 = n�o imprimir)
    integer:: verb_t = 0        !<Indica o "n�vel" de informa��o impressa no terminal durante a otimiza��o

end type HSA_par


type HSA_sai
    logical,private:: inicializado = .false.         !<indica se a vari�vel foi inicializada 
    type(HSA_pto):: ot                      !<Restposta final da otimiza��o
    type(HSA_pto),allocatable:: HM(:)       !<harmony memory
    type(HSA_pto),allocatable:: hAF(:)      !<hist�rico dos pontos em que a fun��o objetivo foi avaliada.
    integer,allocatable:: it_hAF(:)         !<indica a posi��o em que o pto gerado na itera��o est� armazenado em hAF
    integer,allocatable:: ot_hAF(:)         !<indica a posi��o em que o pto �timo na itera��o est� armazenado em hAF
     
    integer:: it = 0        !<itera��o do algoritmo
    integer:: nAF = 0       !<n�mero de avalia��es da fun��o objetivo
    integer:: it_ot = 0     !<N�m. da itera��o em que a solu��o da otimiza��o foi encontrada
    integer:: nAF_ot = 0    !<N�m. de avalia��es da fun��o objetivo at� encontrar a solu��o da otimiza��o
    
    real(8):: tempo         !Tempo total decorrido na otimiza��o (wall clock em segundos)
    
end type HSA_sai


!-----------------------------------------------------------------------------------

    CONTAINS    !-------subrotinas do m�dulo------------

!-----------------------------------------------------------------------------------


subroutine otimiza_HSA(fob,par,sai)
!************************************************************************************
!* PROP�SITO:																		*
!*																					*
!*																					*
!* PROGRAMADOR			: Felipe S. Almeida											*
!* �ltima altera��o em	: agosto/2014												*
!************************************************************************************
!par�metros da subrotina
type(HSA_par),intent(in):: par  !<Par�metros do HSA
type(HSA_sai),intent(out),target:: sai !<Saida do algoritmo

!____________________________________________________________________________________
EXTERNAL fob			! Fun��o externa que calcula o valor da fun��o objetivo
						! CALL FUNC(N_vd,N_vc,X_vd,X_vc,fob)
!____________________________________________________________________________________

    
!vari�veis internas da subrotina                        
integer:: i, j,  iter   

integer:: i_hAF             !<Apontador na hist�ria de AF
integer,pointer:: cont      !<Contador usado para atualizar os par�metros (default sai%it)
integer:: maxcont   !<Contador usado para atualizar os par�metros

integer:: HMS               !<Tamanha atual da mem�ria
integer:: HMvd(par%maxHMS)  !<Componentes da mem�ria para uma vari�vel discreta
real(8):: HMvc(par%maxHMS)  !<Componentes da mem�ria para uma vari�vel cont�nua


real(8):: HMCRvd(maxNvD)    !<HMCR atual para vari�veis discretas
real(8):: PARvd(maxNvD)     !<PAR atual para vari�veis discretas
integer:: PAbw_vd(maxNvD)   !<bw atual para vari�veis discretas

real(8):: HMCRvc(maxNvC)    !<HMCR atual para vari�veis cont�nuas
real(8):: PARvc(maxNvC)     !<PAR atual para vari�veis cont�nuas
real(8):: PAbw_vc(maxNvC)   !<bw atual para vari�veis cont�nuas


logical:: novo_ot           !<Inidicador de obten��o de novo �timo na itera��o

type(HSA_pto):: HI              !<"Harmonia" (novo ponto) criada na itera��o
type(HSA_pto),pointer:: HM(:)   !<Mem�ria em uso
type(HSA_pto),pointer:: ot      !<�timo atual

character(18):: frmt            !formato para escrever dados

real(8):: tempo_inicio
!_____________________________Inicializando Par�metros e vari�veis ________________________

!Obtendo tempo do in�cio da otimiza��o
call cpu_time(tempo_inicio)

!Inicializando vari�veis
call RANDOM_SEED()
call iniciar_HSA_pto(par,HI)
call iniciar_HSA_sai(par,sai)

!definido o contador
if(par%contador=='af'.or.par%contador=='AF') then !N�mero de an�lises da fun��o objetivo como contador
    cont =>sai%nAF
    maxcont = par%maxAF
else !N�mero de itera��es como contador
    cont=>sai%it
    maxcont = par%maxIT
endif
cont = 0

!verifica��o inicial dos par�metros da otimiza��o
if(par%verificar_parametros) call verificar_parametos_HSA_par(par,maxcont)

!incindo o valor dos par�metros de HSA ---
call atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,maxcont,0)
HM=> sai%HM(1:HMS)
ot => sai%ot


!_____________________________Inicializando HM ________________________

sai%HM(:)%fob = huge(0.d0) !Inicia a HM com valor enorme para a fun��o objetivo

if(par%verb_t>0) write(*,'(A)')"HSA>> Forma��o da HM inicial: in�cio"

do i = 1 , HMS !criando os pontos iniciais da HM
    
    do !obtendo um ponto aleatoriamente (verifica se � diferente dos anteriores)
        call rand_HSA_pto(par,HI)
        if(i==1) exit
        if( buscar_HSA_pto(sai%nAF, sai%hAF, HI) == 0 ) exit
    enddo
    
    !avaliando a fun��o objetivo do ponto
    CALL fob(par%nvD, par%nvC, HI%vd, HI%vc, HI%fob )
    
    !armazenando o ponto com fun��o objetivo avaliada
    sai%nAF = sai%nAF+1
    sai%hAF(sai%nAF) = HI
    call atualizar_HM(HI,HMS,HM,j) 
    
    !Definido �timo na itera��o inicial
    if(j==1) then !O novo ponto � o novo �timo
        sai%ot_hAF(1) = sai%nAF
        sai%nAF_ot = sai%nAF
        sai%it_ot = 1
        ot = HI
    endif
enddo

if(par%verb_t>0) write(*,'(A)')"HSA>> Forma��o da HM inicial: fim"


!_____________________________bloco de itera��es do algoritmo________________________

if(par%verb_t>0) write(*,'(A)')"HSA>> Itera��es da otimiza��o: in�cio"

DO !loop das itera��es no algoritmo de busaca pela harmonia 

    
    sai%it = sai%it + 1  !Incrementa o contador de itera��es
    
    
    call atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,maxcont,cont)
    HM=> sai%HM(1:HMS)
   
    !improviso de uma nova harmonia ---- 
    
    !La�o sobre as vari�veis discretas
    IF(par%nvD>0) THEN       
        DO i = 1 , par%nvD
            do j = 1 , HMS ;  HMvd(j) = HM(j)%vd(i) ; enddo !<-transfere os valores da vd(i) de todos os pontos da HM para HMvd(:)
            CALL improvisar_harmonia_vd(HMS, HMCRvd(i), PARvd(i), PAbw_vd(i), (/ par%Li_vD(i), par%Ls_vD(i) /), HMvd, HI%vd(i) )
        END DO 
    END IF
        
    !La�o sobre as vari�veis cont�nuas
    IF(par%nvC>0) THEN       
        DO i = 1 , par%nvC
            do j = 1 , HMS ; HMvc(j) = HM(j)%vc(i) ; enddo !<-transfere os valores da vc(i) de todos os pontos da HM para HMvc(:)
           CALL improvisar_harmonia_vc(HMS, HMCRvc(i), PARvc(i), PAbw_vc(i), (/ par%Li_vC(i), par%Ls_vC(i) /), HMvc, HI%vc(i) )
        END DO 
            
    END IF
        

    !Avalia��o da fun��o objetivo ---
    
    !buscando por pontos avaliados anteriormente (apenas para problemas sem vari�veis cont�nuas)
    i_hAF = 0
    if(par%nvC == 0) i_hAF = buscar_HSA_pto(sai%nAF, sai%hAF, HI)
         
    if(i_hAF/=0) then   !<-para de um ponto que j� foi avaliado anteriormente
        HI%fob = sai%hAF(i_hAF)%fob
        sai%it_hAF(sai%it) = i_hAF
    
    else                !<-para de um ponto que n�o foi avaliado anteriormente
        CALL fob(par%nvD, par%nvC, HI%vd, HI%vc, HI%fob ) !<-avalia��o da fob
        sai%nAF = sai%nAF+1
        sai%hAF(sai%nAF) = HI
        sai%it_hAF(sai%it) = sai%nAF

    endif
    
    
    !Atualia��o da HM ---
    call atualizar_HM(HI,par%maxHMS,sai%HM,j)
    
    
    !Atualiza o valor �timo ---
    novo_ot = j==1 !novo �timo caso HI tenha sido adicionado no in�cio da mem�ria
    
    if(novo_ot) then !Ataliza��o das vari�veis em sai no caso de novo �timo
        ot = HI
        sai%nAF_ot = sai%nAF
        sai%it_ot = sai%it
        sai%ot_hAF(sai%it) = sai%nAF
        
    elseif(sai%it>1) then !Manuten��o de ot_hAF nessa itera��o
        sai%ot_hAF(sai%it) = sai%ot_hAF(sai%it-1)
           
    endif
    
 
    !Imprimindo sa�da da itera��o ---
1001 format("HSA>> cont="I6,"/",I6.6,2x,ES12.4E2)
    if(par%verb_t>0) write(*,1001) cont,maxcont,ot%fob
    
    !Verificar o encerramento das itera��es do HSA  ---
    if(sai%it >= par%maxIT) exit
    if(sai%nAF >= par%maxAF) exit
    !if(cont>=maxcont) exit
    
    
    END DO !Fim do loop de itera��es do HSA   

    call cpu_time(sai%tempo)
    sai%tempo = sai%tempo - tempo_inicio
 
    
    if(par%verb_t>0) then
        write(*,'(A)')"HSA>> Itera��es da otimiza��o: FIM"
        write(*,'(A,es12.4e2)')"HSA>> tempo da otimiza��o: ",sai%tempo
        
        !    write(*,'(A)')"Resultado da otimiza��o:"
    !    if(par%nvd/=0) then
    !        write(*,'(A,i3,A)') "('Vd='," , par%nvd , "(x,i3,','))"
    !        write(frmt,'(A,i3,A)') "('Vd='," , par%nvd , "(x,i3,','))"
    !        write(*,*) frmt
    !        !write(*,frmt)ot%vd(1:par%nvd)
    !    endif
    !    if(par%nvd/=0) then
    !        write(frmt,*)"('Vc=',",par%nvc,"(x,ES12.4E2,','))"
    !        write(*,*) frmt
    !        write(*,frmt)ot%vc(1:par%nvc)
    !    endif
    !    write(frmt,*)"('fob=',x,ES12.4E2)"
    !    write(*,frmt)ot%fob
    endif
    
    call cpu_time(sai%tempo)
    sai%tempo = sai%tempo - tempo_inicio
    
    
    
end subroutine otimiza_HSA
   


 


SUBROUTINE improvisar_harmonia_vd(HMS,HMCR,PAR,PAbw,Lim,HM,x)
!************************************************************************************
!* PROP�SITO:																		*
!* Gerar o valor para uma vari�vel discreta na nova harmonia improvisada            *
!*                                            										*
!************************************************************************************
INTEGER,INTENT(IN):: HMS    !Tamnho da memoria de harmonia (Harmony Memory Size)
REAL(8),INTENT(IN):: HMCR   !Harmony Memory Considering Rate
REAL(8),INTENT(IN):: PAR    !Pitch adjusting rate
INTEGER,INTENT(IN):: PAbw   !Pitch adjusting distance bandwidth
INTEGER,INTENT(IN):: Lim(2) !Limites para o valor gerado (1-m�nimo/2-m�ximo)
INTEGER,INTENT(IN):: HM(HMS)    !Mem�ria de harmonia pra a vari�vel em quest�o
INTEGER,INTENT(OUT):: x     !valor da vari�vel na improvisa��o de um anova harmonia


!vari�veis internas da subrotina
LOGICAL:: teste !auxiliar para decis�es
INTEGER:: ix, ia    !auxiliar

!decidindo se o novo valor ser� gerado usando HM ou ser� gerado aleatoriamente
CALL lRAND(teste, HMCR)

IF(teste) THEN !gerar x usando HM

    !Escolhe aleatoriamente a posi��o em HM par usar na defini��o de x
    CALL iRAND(1,HMS,ix) 
    
    !Define o valor de x
    x = HM(ix)
 
    !decide se ser� realizado o ajuste (pitch adjustment)
    CALL lRAND(teste, PAR)
    
    IF(teste) THEN !realizar of pitch adjustment
        
        !gera o ajuste
        CALL iRAND(1,PAbw,ia)
        
        !decide a dire��o do ajuste na mem�ria
        CALL lRAND(teste) !falso indica  ajuste negativo
        
        !Pitch adjustment  --- levando em conta os limites
        if( x == lim(1) ) then
            x = min(Lim(2), x+ia )
        elseif( x == lim(2) ) then
            x = max(Lim(1), x-ia )
        else
            IF(teste) THEN !ajuste para valor anterior
                x = max(Lim(1), x-ia )
            ELSE
                x = min(Lim(2), x+ia )
            END IF
        endif
              
    END IF          !fim do pitch adjustment

ELSE        !gerar x aleatoriamente

    CALL iRAND(lim(1),lim(2),x)

END IF      !Fim do teste de utiliza��o de HM na gera��o de x     

IA =0    
END SUBROUTINE improvisar_harmonia_vd




SUBROUTINE improvisar_harmonia_vc(HMS,HMCR,PAR,PAbw,Lim,HM,x)
!************************************************************************************
!* PROP�SITO:																		*
!* Gerar o valor para uma vari�vel discreta na nova harmonia improvisada            *
!*                                            										*
!************************************************************************************
INTEGER,INTENT(IN):: HMS    !Tamnho da memoria de harmonia (Harmony Memory Size)
REAL(8),INTENT(IN):: HMCR   !Harmony Memory Considering Rate
REAL(8),INTENT(IN):: PAR    !Pitch adjusting rate
REAL(8),INTENT(IN):: PAbw   !Pitch adjusting distance bandwidth
REAL(8),INTENT(IN):: Lim(2) !Limites para o valor gerado (1-m�nimo/2-m�ximo)
REAL(8),INTENT(IN):: HM(HMS)    !Mem�ria de harmonia pra a vari�vel em quest�o
REAL(8),INTENT(OUT):: x     !valor da vari�vel na improvisa��o de um anova harmonia


!vari�veis internas da subrotina
LOGICAL:: teste !auxiliar para decis�es
INTEGER:: ix    !auxiliar
REAL(8):: ra    !ajuste
REAL(8):: inf, sup  !valores inferior e superior que podem ser usados para PA

!decidindo se o novo valor ser� gerado usando HM ou ser� gerado aleatoriamente
CALL lRAND(teste, HMCR)

IF(teste) THEN !gerar x usando HM

    !Escolhe aleatoriamente a posi��o em HM par usar na defini��o de x
    CALL iRAND(1,HMS,ix) 
    
    !Define o valor de x
    x = HM(ix)
 
    !decide se ser� realizado o ajuste (pitch adjustment)
    CALL lRAND(teste, PAR)
    
    IF(teste) THEN !realizar of pitch adjustment
        
        !gera o ajuste
        CALL rRAND(0.d0,PAbw,ra)
        
        !decide a dire��o do ajuste na mem�ria
        CALL lRAND(teste)
        
        IF(teste) THEN !ajuste para valor anterior
            x = max(Lim(1), x-ra )
        ELSE
            x = min(Lim(2), x+ra )
        END IF
              
    END IF          !fim do pitch adjustment

ELSE        !gerar x aleatoriamente

    CALL rRAND(lim(1),lim(2),x)

END IF      !Fim do teste de utiliza��o de HM na gera��o de x     
    
END SUBROUTINE improvisar_harmonia_vc




subroutine atualizar_HM(pto,HMS,HM,i)
!************************************************************************************
!>Introduz um novo ponto "pto" na Harmony memory "HM" se o valor da fun��o objetivo
!> do novo ponto for 'melhor' (menor valor) do que o pior componente da HM. 
!>A HM fica organizada em ordem crescente do valor da fob.
!>A melhor solu��o atual fica em HM(1) e a pior em HM(HMS)
!************************************************************************************
type(HSA_pto),intent(in)::    pto       !<Novo vetor solu��o
integer,intent(in):: HMS                !<Tamanho de HM
type(HSA_pto),intent(inout):: HM(:)     !<Harmony memory
integer,intent(out)::        i          !<Posi��o do novo vetor na HM (obs: i=0 se pto n�o for adicionado)

integer:: j, k
i=0 !Inicia com posi��o "0" indicando que a HM n�o foi atualizada.
if(pto%fob>HM(HMS)%fob) return !Sai caso pto seja pior que todos na HM
do j = 1 , HMS 
    if(pto%fob<HM(j)%fob) then !Encontrada a posi��o para inserir pto
        if(j/=HMS) then        !Deslocando pontos piores em uma posi��o
            do k = HMS , j+1 , -1 
                HM(k) = HM(k-1)
            enddo
        endif
        i=j                     !Posi��o de pto na HM
        HM(i)=pto               !Adicionando pto na HM
        exit                    !sai do la�o quando encontrou a posi��o
    endif
enddo

end subroutine



!---------------------- rotinas relacionadas � vari�vel tipo HSA_par ----------------



subroutine atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,imax,ii)
!************************************************************************************
!>Retorna o valor das vari�veis com base no valor do contador atual
!************************************************************************************
type(HSA_par),intent(in):: par              !<Par�metros informados inicialmente
integer,intent(out):: HMS                   !Tamanho da mem�ria usada (HMS)
real(8),intent(out):: HMCRvd(:), PARvd(:)   !<HMCR e PAR para vari�veis discretas
integer,intent(out):: PAbw_vd(:)            !bw para vari�veis discretas
real(8),intent(out):: HMCRvc(:), PARvc(:)   !<HMCR e PAR par vari�veis cont�nuas
real(8),intent(out):: PAbw_vc(:)            !<bw para vari�veis cont�nuas
integer,intent(in):: imax                   !Valor limite do contador
integer,intent(in):: ii                     !<Contador

integer:: jj, tt
real(8):: rr, mm        !valor real de ii e imax, respectivamente
character(30) err_str !string com mensagem de erro

!vari�veis para controle da impress�o dos par�metos em arquivo
logical:: arq_aberto !informa��o se o arquivo est� aberto
character(50):: arq_action, frmt1

!----
300 format("ERR> @atualizar_parametros_HSA_par: em ii=",i9,A30)

    
err_str = ' '
rr = real(ii,8) !converte ii em real para ser usado na fun��o fnc
mm = real(imax,8)

!---

HMS = nint( fnc(par%HMS_tf,par%HMS_c,mm,rr) )
if(HMS>par%maxHMS) then 
    write(err_str,"(' HMS = ',i3,' > maxHMS = ',i3)")HMS,par%maxHMS
    goto 301
endif

!---

if(par%nvD>0) then !par�metros para as vari�veis discretas
    do jj = 1 ,  par%nvD
        HMCRvd(jj) = fnc( par%HMCRvd_tf(jj) , par%HMCRvd_c(:,jj) ,mm ,rr )
        PARvd(jj)  = fnc( par%PARvd_tf(jj)  , par%PARvd_c(:,jj)  ,mm ,rr )
        PAbw_vd(jj) = nint( fnc( par%PAbw_vd_tf(jj) , par%PAbw_vd_c(:,jj) ,mm, rr ) )
 
        if( HMCRvd(jj)<0.d0 .or. 1.d0<HMCRvd(jj))then
            write(err_str,"(' HMCRvd(',i3,')= ',f6.2)")jj,HMCRvd(jj); goto 301
        else if( PARvd(jj)<0.d0 .or. 1.d0<PARvd(jj)) then
            write(err_str,"('  PARvd(',i3,')= ',f6.2)")jj,PARvd(jj); goto 301
        endif
        
    enddo
endif

!---

if(par%nvC>0) then !par�metros para as vari�veis cont�nuas
    do jj = 1 ,  par%nvC
        HMCRvc(jj) = fnc( par%HMCRvc_tf(jj) , par%HMCRvc_c(:,jj) ,mm ,rr )
        PARvc(jj)  = fnc( par%PARvc_tf(jj)  , par%PARvc_c(:,jj)  ,mm ,rr )
        PAbw_vc(jj) = fnc( par%PAbw_vc_tf(jj) , par%PAbw_vc_c(:,jj) ,mm ,rr )

        if( HMCRvc(jj)<0.d0 .or. 1.d0<HMCRvc(jj))then
            write(err_str,"(' HMCRvc(',i3,')= ',f6.2)")jj,HMCRvc(jj); goto 301
        else if( PARvc(jj)<0.d0 .or. 1.d0<PARvc(jj)) then
            write(err_str,"('  PARvc(',i3,')= ',f6.2)")jj,PARvc(jj); goto 301
        endif
        
    enddo
endif

!---

!!impress�o do valor dos par�metros em arquivo
call imprimir_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,ii,par%uni_imp_par)

!---

return !retorna sem problemas

301 write(*,300)ii,err_str !informa sobre o problema encontrado

end subroutine atualizar_parametros_HSA_par



subroutine imprimir_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,ii,uni)
!************************************************************************************
!>Imprime o valor das vari�veis em um arquivo
!************************************************************************************
type(HSA_par),intent(in):: par              !<Par�metros informados inicialmente
integer,intent(in):: HMS                   !Tamanho da mem�ria usada (HMS)
real(8),intent(in):: HMCRvd(:), PARvd(:)   !<HMCR e PAR para vari�veis discretas
integer,intent(in):: PAbw_vd(:)            !bw para vari�veis discretas
real(8),intent(in):: HMCRvc(:), PARvc(:)   !<HMCR e PAR par vari�veis cont�nuas
real(8),intent(in):: PAbw_vc(:)            !<bw para vari�veis cont�nuas
integer,intent(in):: ii                     !<Contador
integer,intent(in):: uni                    !Unidade para a impress�o dos par�metros


integer:: jj, tt
character(30) err_str !string com mensagem de erro

!vari�veis para controle da impress�o dos par�metos em arquivo
logical:: arq_aberto !informa��o se o arquivo est� aberto
character(50):: arq_action, frmt1

!----
300 format("ERR> @imprimir_parametros_HSA_par: em ii=",i9,A30)

err_str = ' '

!impress�o do valor dos par�metros em arquivo
if(uni/=0) then

    if(ii==0) then !a��es no in�cio da otimiza��o
        
        !verificando a conformidade do arquivo associado � unidade uni
        inquire(unit=uni, opened= arq_aberto,action=arq_action)
        if(.not. arq_aberto) then
            write(err_str,'(A,i5)') "Arquivo n�o aberto, unidade=" , uni
            goto 301
        else if(arq_action=='READ') then
            write(err_str,'(A,i5)') "Arquivo somente para leitura, unidade=" , uni
            goto 301
        endif        
        
        !cabe�alho com lista de vari�veis (formato tecplot)  
        write(uni,'(A)',advance='NO') "VARIABLES = it, HMS"
        if(par%NvD>0)then
            do jj = 1 , par%NvD
                write(uni,'(", HMCR-vD",i3.3)',advance="NO")jj
            enddo
            do jj = 1 , par%NvD
                write(uni,'(", PAR-vD",i3.3)',advance="NO")jj
            enddo
            do jj = 1 , par%NvD
                write(uni,'(", bw-vD",i3.3)',advance="NO")jj
            enddo
        endif
        if(par%NvC>0)then
            do jj = 1 , par%NvC
                write(uni,'(", HMCR-vC",i3.3)',advance="NO")jj
            enddo
            do jj = 1 , par%NvC
                write(uni,'(", PAR-vC",i3.3)',advance="NO")jj
            enddo
            do jj = 1 , par%NvC
                write(uni,'(", bw-vC",i3.3)',advance="NO")jj
            enddo
        endif
         write(uni,*)""
    endif
    
    !imprimindo par�metros
    jj = par%NvD
    tt = par%NvC
    if(par%NvD>0 .and. par%NvC>0) then
        write(frmt1,'(A,4(i3,A))') '(i6,x,i3,',2*jj,'(x,f5.3),',jj,'(x,i3),',2*tt,'(x,f5.3),',tt,'(x,ES12.3))'
        write(uni,frmt1) ii, HMS, HMCRvd(1:jj), PARvd(1:jj), PAbw_vd(1:jj), &
                                              HMCRvc(1:tt), PARvc(1:tt), PAbw_vc(1:tt)
    else if(par%NvD>0) then
        write(frmt1,'(A,2(i3,A))') '(i6,x,i3,',2*jj,'(x,f5.3),',jj,'(x,i3))'
        write(uni,frmt1) ii, HMS, HMCRvd(1:jj), PARvd(1:jj), PAbw_vd(1:jj)
    else
        write(frmt1,'(A,2(i3,A))') '(i6,x,i3,',2*tt,'(x,f5.3),',tt,'(x,ES12.3))'
        write(uni,frmt1) ii, HMS, HMCRvc(1:jj), PARvc(1:jj), PAbw_vc(1:jj)
    endif
     
endif

return

301 write(*,300)ii,err_str !informa sobre o problema encontrado

end subroutine imprimir_parametros_HSA_par


subroutine verificar_parametos_HSA_par(par,imax)
!************************************************************************************
!>	Verifica se haver� problemas com os par�metros gerados durante a otimiza��o
!************************************************************************************
type(HSA_par),intent(in):: par              !<Par�metros informados inicialmente
integer,intent(in):: imax                   !Valor limite do contador


integer:: HMS                       !Tamanho da mem�ria usada (HMS)
real(8):: HMCRvd(maxNvD), PARvd(maxNvD)   !<HMCR e PAR para vari�veis discretas
integer:: PAbw_vd(maxNvD)              !bw para vari�veis discretas
real(8):: HMCRvc(maxNvC), PARvc(maxNvC)   !<HMCR e PAR par vari�veis cont�nuas
real(8):: PAbw_vc(maxNvC)              !<bw para vari�veis cont�nuas
integer:: ii                        !<Contador

character(1):: resp
integer:: uni , err
logical:: aberto

!Abrindo um arquivo para imprimir a verifica��o dos dados
uni=0
do
    uni = uni + 1
    inquire(unit=uni,opened=aberto)
    if(.not. aberto) open(unit=uni,file="verificar_parametos_HSA_par.plt",action="WRITE",iostat=err)
    if(err==0) exit
    if(err==500) then
        write(*,*)"ERR>> @verificar_parametos_HSA_par : n�o foi poss�vel abrir um arquivo para imprimir &
        o arquivo de verifica��o"
        goto 401
    endif
    
enddo

if(par%verb_t>0) write(*,*)"@verificar_parametos_HSA_par : Verificando os par�metros da otimiza��o"

!Gerando e verificando os valores dos par�metros que ser�o usados durante a otimiza��o
do ii = 0 , imax
    call atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,imax,ii)
    call imprimir_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,ii,uni)
enddo

401 write(*,'(A)',advance='NO')"Prosseguir com a otimiza��o? (S/N):"
read(*,*)resp
if(resp=='S' .or. resp=='s') then
    write(*,*) "Iniciando a otimiza��o ... "
else
    write(*,*) "Abortando a otimiza��o ... "
    stop
endif

close(uni)

end subroutine verificar_parametos_HSA_par



real(8) function fnc(tf,c,m,x)
!************************************************************************************
!>	Retorna o valor de uma das fun��es programadas 
!> O tipo de fun��o � defindi por tf.
!> O valor � retornado para vari�vel independente x, que deve estar entre 0 e m
!> m � o valor m�ximo de x
!************************************************************************************
integer,parameter:: nnc = ncoef     !<n�mero de constantes
integer,parameter:: nntf = nctf     !<n�mero de caracteres definindo a fun��o
!Argumentos da fun��o========================================
character(nntf),intent(in)::   tf   !<tipo de fun��o
real(8),intent(in)::        c(nnc)  !<Coeficientes da fun��o
real(8),intent(in)::        m       !<Valor m�ximo de x
real(8),intent(in)::        x       !<Vari�vel independente

selectcase(tf) !selecionando o tipo de fun��o
case('cst') !constante
    fnc = c(1)
case('lin') !linear ; fnc(0)=C1 e fnc(m)=C2
    fnc = c(1) + (c(2)-c(1))*x/m
case('pvf') !par�bola com v�rtice no final; fnc(0)=c1 ; fnc(m)=c2 ; df(m)/dx = 0
    fnc = c(1) + (c(2)-c(1))* (x/m) * (2 - x/m) 
case('pv0') !par�bola com v�rtice no in�cio; fnc(0)=c1 ; fnc(m)=c2 ; df(0)/dx = 0
    fnc = c(1) + ( c(2)-c(1))*(x/m)**2
case('pim') !par�bola com in�cio e m�ximo/m�nimo definido; fnc(0)=c1; fnc(c2)=c3; df(c2)/dx=0
    fnc = c(1) + (c(3)-c(1))* (x/c(2)) * (2 - x/c(2))
case('pg1') !polin�mio do primeiro grau (linear)
    fnc = c(1) + c(2)*x
case('pg2') !polin�mio do segundo grau (par�bola)
    fnc = c(1) + c(2)*x + c(3)*x**2
case('pg3') !polin�mio do terceiro grau (c�bica)
    fnc = c(1) + c(2)*x + c(3)*x**2 + c(4)*x**3
case('exp') !exponencial
    fnc = c(1) * exp( c(2) * x ) + c(3)
case default    
    fnc = huge(fnc)
    write(*,1001) tf
1001 format("ERR>> @fnc: tipo de funcao nao encontrada = ",A3)     
endselect

end function fnc



!---------------------- rotinas relacionadas � vari�vel tipo HSA_pto ----------------

subroutine iniciar_HSA_pto(par,pto)
!************************************************************************************
!>Aloca as vari�veis vd e vc de pto com tamanhos definidos em par%nvd e par%nvc
!>, al�m de zerar dotas as sub-vari�veis de pto
!************************************************************************************
type(HSA_par),intent(in):: par  !<Par�metros
type(HSA_pto),intent(out):: pto !<ponto no espa�o de resposta

integer:: sa, err

if( allocated(pto%vd) ) deallocate(pto%vd)
if(par%nvD>0) then
    allocate(pto%vd(par%nvD),stat=sa)
    pto%vd=0
endif

if( allocated(pto%vc) ) deallocate(pto%vc)
if(par%nvC>0) then
    allocate(pto%vc(par%nvC),stat=sa)
    pto%vc=0.d0
endif

pto%fob = 0.d0

end subroutine iniciar_HSA_pto



subroutine terminar_HSA_pto(pto)
!************************************************************************************
!> Desaloca as vari�veis vd e vc e zera o valor de fob
!************************************************************************************
type(HSA_pto),intent(inout):: pto !<ponto no espa�o de resposta
if(allocated(pto%vd)) deallocate(pto%vd)
if(allocated(pto%vc)) deallocate(pto%vc)
pto%fob = 0.d0
end subroutine terminar_HSA_pto



    
subroutine rand_HSA_pto(par,pto)
!************************************************************************************
!>Gera valores aleat�rios para as vari�ves vd e vc, dentro dos limites estabelecidos
!> na vari�vel par%Li_** e par%Li_** correspondente.
!************************************************************************************
type(HSA_par),intent(in):: par      !<Par�metros
type(HSA_pto),intent(inout):: pto   !<ponto no espa�o de resposta

if(par%nvD>0) call iRAND( par%Li_vD , par%Ls_vD, pto%vd) 
if(par%nvC>0) call rRAND( par%Li_vC , par%Ls_vC, pto%vc) 

end subroutine rand_HSA_pto    
    


function comparar_HSA_pto(pr,pc) result(comp)
!************************************************************************************
!> Compara dois HSA_pto, pc e pr, apresentando resutados em comp(i),
!> i=1 p/ var. discretas e i=2 para vari�veis cont�nuas            
!*																					
!>	comp(i)= 0 - indica igualdade													
!>	comp(i)< 0 - indica que pc tem a vari�veil de posi��o abs(comp(i)) menor que em pr		
!>	comp(i)> 0 - indica que pc tem a vari�veil de posi��o comp(i) maior que em pr		
!>  comp(i)= huge - indica que pr e pc n�o tem o mesmo n�mero de vari�veis desse tipo
!> 	
!************************************************************************************
!par�metros da subrotina
type(HSA_pto),intent(in):: pr   !<ponto de refer�ncia
type(HSA_pto),intent(in):: pc   !<ponto comparado com a refer�ncia
integer,dimension(2):: comp     !<vetor indicando a posi��o em que existe diferen�a

integer:: i, n!vari�veis internas

comp = 0

!--- vari�veis discretas
!N�mero de vari�veis discertas 
n = size(pr%vd)
if(n/=size(pc%vd)) then !No caso de ter 
    comp(1) = huge(i)
else if(n>0) then      !No caso de haver vari�veis discretas
    do i = 1 , n
        if(pc%vd(i)==pr%vd(i)) cycle
        comp(1) = sign(i, pc%vd(i)-pr%vd(i) )
        exit
    enddo
endif


!--- vari�veis cont�nuas
!N�mero de vari�veis cont�nuas em cada ponto
n = size(pr%vc)
if(n/=size(pc%vc)) then !No caso de ter 
    comp(2) = huge(i)
else if(n>0) then      !No caso de haver vari�veis discretas
    do i = 1 , n
        if(pc%vc(i)==pr%vc(i)) cycle
        comp(2) = i
        if(pc%vc(i)<pr%vc(i)) comp(2) = -i
        exit
    enddo
endif
    
end function comparar_HSA_pto




function buscar_HSA_pto(nAF,hAF,pto) result(pos)
!************************************************************************************
!< Busca o ponto pto dentre os pontos armazenados em hAF
!< Retorna: 
!< pos = 0 se pto n�o for encontrado em hAF 
!< pos igual � posi��o de pto em hAF							
!************************************************************************************
!par�metros da subrotina
integer,intent(in):: nAF            !<n�mero de pontos armazenados em hAF
type(HSA_pto),intent(in):: hAF(:)   !<pontos analisados anteriormente
type(HSA_pto),intent(in):: pto      !<ponto buscado

integer:: pos !<posi��o de pto em hAF
integer:: comp(2)
integer:: i

pos = 0

do i = 1 , nAF
    comp = comparar_HSA_pto(hAF(i),pto)
    if(comp(1)/=0 .or. comp(2) /= 0 ) cycle
    pos=i
    exit
enddo
  
end function buscar_HSA_pto




subroutine escrever_formatado_HSA_pto(dtv, unit, iotype, v_list, iostat, iomsg)
!************************************************************************************
!>Impress�o formatada da vari�vel HSA_pto
!************************************************************************************
type(HSA_pto), intent(in) :: dtv
integer, intent(in) :: unit
character(len=*), intent(in) :: iotype
integer, intent(in) :: v_list(:)
integer, intent(inout) :: iostat
integer, intent(inout) :: iomsg

integer:: sl !n�mero de informa��es passadas em v_list
integer:: nvd, lvd !n�mero de vari�veis discretas e n�m. de d�gitos para impress�o
integer::nvc !N�mero de vari�veis cont�nuas
integer:: f  !C�digo para o formato de impress�o
character(50):: frmt !formato para impress�o
integer:: ii    !vari�veis auxiliares
 
sl = size(v_list)
lvd = 4 !Adotando vari�veis inteiras com at� 4 d�gitos

selectcase(sl)
case(:-1) ; iostat = 1 !; iomsg = 'v_list'!'size(v_list)<0'
case(0);  f = 1 ; nvd = size(dtv%vd) ; nvc = size(dtv%vc)
case(1:2);  f = v_list(1) ; nvd = size(dtv%vd) ; nvc = size(dtv%vc)
case(3);   f = v_list(1) ; nvd = v_list(2) ; nvc = v_list(3)
case(4:);  f = v_list(1) ; nvd = v_list(2) ; nvc = v_list(3) ; lvd = v_list(4)
case default ; f = 1 ; nvd = size(dtv%vd) ; nvc = size(dtv%vc)
end select

        
selectcase(f)
case(1)
    write(frmt,*)  "(ES13.6,x,",nvd,"(i",lvd,",',')",nvc,"(x,ES13.6,','))"
case(2)
    write(frmt,*)  "(fob = ES13.6,/,'vD = '",nvd,"(i",lvd,",','),/,' vC = '",nvc,"(x,ES13.6,','))"
case default
    
    write(frmt,*)  "(ES13.6,x," , nvd , "(i" , lvd , ",',')" , nvc , "(x,ES13.6,','))"
    
end select
    
    write(unit,frmt) dtv%fob, dtv%vd, dtv%vc
    
end subroutine escrever_formatado_HSA_pto
    
    


!---------------------- rotinas relacionadas � vari�vel tipo HSA_sai ----------------

    
subroutine iniciar_HSA_sai(par,sai)
!************************************************************************************
!* PROP�SITO:																		*
!*																					*
!*																					*
!* PROGRAMADOR			: Felipe S. Almeida											*
!* �ltima altera��o em	: agosto/2014												*
!************************************************************************************
type(HSA_par):: par
type(HSA_sai):: sai

integer:: i

if(sai%inicializado) call terminar_HSA_sai(sai)

sai%inicializado = .true.

call iniciar_HSA_pto(par,sai%ot)

allocate( sai%HM(par%maxHMS) )
do i = 1 , par%maxHMS
    call iniciar_HSA_pto(par,sai%HM(i))
enddo

allocate( sai%hAF(par%maxAF) )
do i = 1 , par%maxAF
    call iniciar_HSA_pto(par,sai%hAF(i))
enddo

allocate( sai%it_hAF(par%maxIT) ) 
allocate( sai%ot_hAF(par%maxIT) ) 


sai%it = 0
sai%nAF = 0
sai%it_ot = 0
sai%nAF_ot = 0
sai%tempo = 0.d0

end subroutine iniciar_HSA_sai

    
subroutine terminar_HSA_sai(sai)
!************************************************************************************
!>Desaloca as vari�veis alocaveis de HSA_sai e zera as demais vari�veis
!************************************************************************************
type(HSA_sai),intent(inout):: sai

integer:: i

sai%inicializado = .false.         

call terminar_HSA_pto(sai%ot)    

do i = 1 , size(sai%HM)
    call terminar_HSA_pto(sai%HM(i))
enddo

do i = 1 , size(sai%hAF)
    call terminar_HSA_pto(sai%hAF(i))
enddo

deallocate(sai%HM)
deallocate(sai%hAF)
deallocate(sai%it_hAF)
deallocate(sai%ot_hAF)

sai%it = 0
sai%nAF = 0
sai%it_ot = 0
sai%nAF_ot = 0
sai%tempo = 0.d0

end subroutine terminar_HSA_sai




!subroutine imprimir_HSA_sai(par,sai,uni)
!!************************************************************************************
!!>Imprime as informa��es sobre o processo de otimiza��o em um arquivo
!!************************************************************************************
!type(HSA_par),intent(in):: par  !Par�metros da otimiza��o
!type(HSA_sai),intent(in)::  sai !Dados resultantes da otimiza��o
!integer,intent(in):: uni        !Unidade para impress�o dos dados contidos em sai
!
!!vari�veis internas
!character(len=15):: rsp
!logical:: aberto
!integer:: ii, lvd
!character(300):: txt
!
!!----Verificando o arquivo atrelato � unidade uni
!inquire(unit=uni,opened=aberto,write=rsp)
!if(.not.aberto .or. rsp=="NO") then
!    write(*,*)"ERR>>@imprimir_HSA_sai: n�o � poss�vel imprimir no arquivo"
!endif
!
!----
!ii = max( maxval( abs(par%Li_vd) ) , maxval( abs(par%Ls_vd)) )
!lvd = ceiling( log10( max( 9.d0 ,  real(ii , 8 ) ) ) ) ; pause
!
!end subroutine imprimir_HSA_sai




subroutine imprimir_HSA_sai_hAF_plt(par,sai,uni)
!************************************************************************************
!>Imprime as informa��es sobre o processo de otimiza��o em um arquivo
!************************************************************************************
type(HSA_par),intent(in):: par  !Par�metros da otimiza��o
type(HSA_sai),intent(in)::  sai !Dados resultantes da otimiza��o
integer,intent(in):: uni        !Unidade para impress�o dos dados contidos em sai

!vari�veis internas
integer:: jj, i, nit, j1,j2
character(300):: txt
character(100):: frmt

integer,allocatable:: aVd(:)
real(8),allocatable:: aVc(:)
real(8):: afob


write(uni,'(2/,A)') 'TITLE = "Historico HSA"'

!Linha com o nome das vari�veis:
write(uni,'(2/,A)',advance='NO')'VARIABLES = E, fob,'
if(par%NvD>0)then
    do jj = 1 , par%NvD
        write(uni,'(", Vd",i3.3)',advance="NO")jj
    enddo
endif
if(par%NvC>0)then
    do jj = 1 , par%NvC
        write(uni,'(", Vc",i3.3)',advance="NO")jj
    enddo
endif       
write(frmt,*) '(i5,2x,ES15.8,2x,',par%NvD,'(i3,2x),',par%NvC,'(ES15.8,2x))'


nit = sai%it_ot

!write(uni,'(2/,A,i7,a)') 'ZONE T="it_aAF", I=',nit,', F=POINT'
!do i = 1 , nit
!    jj = sai%it_hAF(i)
!    write(uni,frmt)i,sai%hAF(jj)%fob,sai%hAF(jj)%Vd(:),sai%hAF(jj)%Vc(:)
!enddo


!write(uni,'(2/,A,i7,a)') 'ZONE T="ot_hAF", I=',nit,', F=POINT'
write(uni,'(2/,A,i7,a)') 'ZONE T="ot_hAF", F=POINT'
j1=0
j2=0
do i = 1 , nit
    jj = sai%ot_hAF(i)
    if(j1==jj)cycle
    j1=jj
    j2=j2+1
    write(uni,frmt)j2,sai%hAF(jj)%fob,sai%hAF(jj)%Vd(:),sai%hAF(jj)%Vc(:)
enddo


!allocate(aVd(par%nVd),aVc(par%nVc))
!
!write(uni,'(2/,A,i7,a)') 'ZONE T="ot_hAF", I=',nit,', F=POINT'
!aVd  = 0
!aVc  = 0.d0
!afob = 0.d0
!write(uni,frmt)i,afob,avd,avc
!
!do i = 2 , nit
!    j1 = sai%ot_hAF(i)
!    j2 = sai%ot_hAF(i-1)
!    
!    aVd  = aVd  + abs( sai%hAF(j1)%Vd(:) - sai%hAF(j2)%Vd(:) )
!    aVc  = aVc  + abs( sai%hAF(j1)%Vc(:) - sai%hAF(j2)%Vc(:) )
!    afob = afob + abs( sai%hAF(j1)%fob - sai%hAF(j2)%fob )
!    write(uni,frmt)i,afob,avd,avc
!enddo



end subroutine imprimir_HSA_sai_hAF_plt



END MODULE HSAg