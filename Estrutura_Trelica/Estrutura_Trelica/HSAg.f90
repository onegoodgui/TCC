MODULE HSAg
!************************************************************************************
!* PROPÓSITO:																		*
!* Realiza otimização pelo Algoritmo de Busca pela Harmonia (HSA)					*
!* (HSA = Harmony Search algorithm)           										*
!* 																					*
!* Solução de problemas de otimização (lineares ou não-lineares)					*
!* Variáveis inteiras ou variáveis discretas conectadas a inteiras.					*
!* (linked-discrete variables nonlinear programing LD-NLP)							*
!* Variáveis reais																	*
!* 																					*
!* 																					*
!*	O problema de otimização deve ser do tipo:										*
!*																					*
!*		Minimizar: F(X)																*
!*		tal que X é um vetor com n compnentes x_i									*
!* 		e x_i_L <= x_i <= x_i_U														*
!* 																					*
!*      F(X) é a função objetivo													*
!* 		No caso de otimizações com restrições F(X) deve contemplar funções de 		*
!* 		penalidade																	*
!* 																					*
!* Referências:  																	*
!*		1) Z.W.Geem,J.-H.Kim,G.V.Loganathan,A new heuristic optimization            *
!*		algorithm: harmony search.Simulation 76(2)(2001)60–68.	                    *
!* 		2) A new meta-heuristic algorithm for continuous engineering optimization   *
!* 		 - harmony search theory and practice.                                      *
!* 		 Comput.Methods Appl. Mech. Engrg.194(2005)3902–3933						*
!* 		3) Optimum design of steel frames using a harmony search algorithm.         * 
!*      Struct. Multidisc Optim (2008) 36:393-401                                   *
!* 																					*
!* 																					*
!* 																					*
!*																					*
!* Autor: Felipe Schaedler de Almeida (Almeida, F. S.)	 			
!************************************************************************************
USE numRAND     !Módulo numRAND com subrotinas auxiliares para geração de valores
                !randômicos lógicos, inteiros ou reais em intervalos definidos

IMPLICIT NONE   !Todas as variáveis devem ser declaradas
PRIVATE         !Todas as rotainas são privadas a menos da lista abaixo

!-------------Lista de subrotinas e variáveis públicas-------------
public otimiza_HSA
public iniciar_HSA_sai
public terminar_HSA_sai

public imprimir_HSA_sai_hAF_plt


public HSA_pto
public HSA_par
public HSA_sai


!public fnc


!----------Parâmetros internos do algoritmo  -----------------
!(usadas para alocação estática de algumas variáveis)
integer,parameter:: maxNvD = 50 !<Núm. máximo de variáveis discretas
integer,parameter:: maxNvC = 50 !<Núm. máximo de variáveis contínuas
integer,parameter:: ncoef = 5   !<Núm. de coeficientes das funções de variação do parâmetros
integer,parameter:: nctf = 3    !<Núm. de caracteres usados para definir as funções 


!-----------Declaração das variáveis tipo derivadas--------------

type HSA_pto !<ponto avaliado pelo algoritmo no espaço de respostas do problema
    integer,allocatable:: vd(:) !<variável discreta
    real(8),allocatable:: vc(:) !<variáveis contínuas
    real(8):: fob = 0.d0        !<Valor da função objetivo
end type HSA_pto

type HSA_par !<Parâmetros do HSA
    
    integer:: nvD  = 0          !<número de variáveis discretas
    integer:: nvC  = 0          !<número de variáveis contínuas
    integer:: Li_vd(maxNvD) = 0 !<Limite inferior das variáveis discretas
    integer:: Ls_vd(maxNvD) = 0 !<Limite superior das variáveis discretas
    real(8):: Li_vc(maxNvC) = 0.d0    !<Limite inferior das variáveis contínuas
    real(8):: Ls_vc(maxNvC) = 0.d0    !<Limite superior das variáveis contínuas  
    
    integer:: maxIT = huge(1)   !Núm. máximo de iterações do algoritmo
    integer:: maxAF = huge(1)   !Núm. máximo de avaliações da função objetivo (sem considerar as repetições)
    
    logical:: verificar_parametros = .false. !Indica se deve ser realizada a verificação dos parâmetros previamente à otimização

    character(2):: contador = 'it'      !<define se usa o núm. de iterações ('it' e default) ou o núm de análises da função objetivo (af) como contador
    integer:: maxHMS = 0                !<Tamanho máximo da HM
    character(nctf):: HMS_tf  = 'cst'   !<Função de variação de HM
    real(8):: HMS_c(ncoef) = 0          !<Coeficientes da função de variação de HM
    
    !Parâmetros para as variáveis discretas  ---vd---    
    character(nctf):: HMCRvd_tf(maxNvD) = 'cst'   !<Função de variação de HMCR para vd
    real(8):: HMCRvd_c(ncoef,maxNvD) = 0.d0       !<Coef. da função de variação de HMCR para vd

    character(nctf):: PARvd_tf(maxNvD) = 'cst'    !<Função de variação de PAR para vd
    real(8):: PARvd_c(ncoef,maxNvD) = 0.d0        !<Coef. da função de variação de PAR para vd
    
    character(nctf):: PAbw_vd_tf(maxNvD) = 'cst'  !<Função de variação de bw para vd
    real(8):: PAbw_vd_c(ncoef,maxNvD) = 0         !<Coef. da função de variação de bw para vd
    
    !Parâmetros para as variáveis constínuas  ---vc---    
     
    character(nctf):: HMCRvc_tf(maxNvC) = 'cst'   !<Função de variação de HMCR para vc
    real(8):: HMCRvc_c(ncoef,maxNvC) = 0.d0       !<Coef. da função de variação de HMCR para vc
    
    character(nctf):: PARvc_tf(maxNvC) = 'cst'    !<Função de variação de PAR para vc
    real(8):: PARvc_c(ncoef,maxNvC) = 0.d0        !<Coef. da função de variação de PAR para vc
    
    character(nctf):: PAbw_vc_tf(maxNvC) = 'cst'  !<Função de variação de bw para vc
    real(8):: PAbw_vc_c(ncoef,maxNvC) = 0.d0      !<Coef. da função de variação de bw para vc
    
    !Parâmetros para controle de impressão de dados durante a otimização
    integer:: uni_imp_par = 0   !<Unid. para impressão de parâmetros do HSA (0 = não imprimir)
    integer:: verb_t = 0        !<Indica o "nível" de informação impressa no terminal durante a otimização

end type HSA_par


type HSA_sai
    logical,private:: inicializado = .false.         !<indica se a variável foi inicializada 
    type(HSA_pto):: ot                      !<Restposta final da otimização
    type(HSA_pto),allocatable:: HM(:)       !<harmony memory
    type(HSA_pto),allocatable:: hAF(:)      !<histórico dos pontos em que a função objetivo foi avaliada.
    integer,allocatable:: it_hAF(:)         !<indica a posição em que o pto gerado na iteração está armazenado em hAF
    integer,allocatable:: ot_hAF(:)         !<indica a posição em que o pto ótimo na iteração está armazenado em hAF
     
    integer:: it = 0        !<iteração do algoritmo
    integer:: nAF = 0       !<número de avaliações da função objetivo
    integer:: it_ot = 0     !<Núm. da iteração em que a solução da otimização foi encontrada
    integer:: nAF_ot = 0    !<Núm. de avaliações da função objetivo até encontrar a solução da otimização
    
    real(8):: tempo         !Tempo total decorrido na otimização (wall clock em segundos)
    
end type HSA_sai


!-----------------------------------------------------------------------------------

    CONTAINS    !-------subrotinas do módulo------------

!-----------------------------------------------------------------------------------


subroutine otimiza_HSA(fob,par,sai)
!************************************************************************************
!* PROPÓSITO:																		*
!*																					*
!*																					*
!* PROGRAMADOR			: Felipe S. Almeida											*
!* Última alteração em	: agosto/2014												*
!************************************************************************************
!parâmetros da subrotina
type(HSA_par),intent(in):: par  !<Parâmetros do HSA
type(HSA_sai),intent(out),target:: sai !<Saida do algoritmo

!____________________________________________________________________________________
EXTERNAL fob			! Função externa que calcula o valor da função objetivo
						! CALL FUNC(N_vd,N_vc,X_vd,X_vc,fob)
!____________________________________________________________________________________

    
!variáveis internas da subrotina                        
integer:: i, j,  iter   

integer:: i_hAF             !<Apontador na história de AF
integer,pointer:: cont      !<Contador usado para atualizar os parâmetros (default sai%it)
integer:: maxcont   !<Contador usado para atualizar os parâmetros

integer:: HMS               !<Tamanha atual da memória
integer:: HMvd(par%maxHMS)  !<Componentes da memória para uma variável discreta
real(8):: HMvc(par%maxHMS)  !<Componentes da memória para uma variável contínua


real(8):: HMCRvd(maxNvD)    !<HMCR atual para variáveis discretas
real(8):: PARvd(maxNvD)     !<PAR atual para variáveis discretas
integer:: PAbw_vd(maxNvD)   !<bw atual para variáveis discretas

real(8):: HMCRvc(maxNvC)    !<HMCR atual para variáveis contínuas
real(8):: PARvc(maxNvC)     !<PAR atual para variáveis contínuas
real(8):: PAbw_vc(maxNvC)   !<bw atual para variáveis contínuas


logical:: novo_ot           !<Inidicador de obtenção de novo ótimo na iteração

type(HSA_pto):: HI              !<"Harmonia" (novo ponto) criada na iteração
type(HSA_pto),pointer:: HM(:)   !<Memória em uso
type(HSA_pto),pointer:: ot      !<Ótimo atual

character(18):: frmt            !formato para escrever dados

real(8):: tempo_inicio
!_____________________________Inicializando Parâmetros e variáveis ________________________

!Obtendo tempo do início da otimização
call cpu_time(tempo_inicio)

!Inicializando variáveis
call RANDOM_SEED()
call iniciar_HSA_pto(par,HI)
call iniciar_HSA_sai(par,sai)

!definido o contador
if(par%contador=='af'.or.par%contador=='AF') then !Número de análises da função objetivo como contador
    cont =>sai%nAF
    maxcont = par%maxAF
else !Número de iterações como contador
    cont=>sai%it
    maxcont = par%maxIT
endif
cont = 0

!verificação inicial dos parâmetros da otimização
if(par%verificar_parametros) call verificar_parametos_HSA_par(par,maxcont)

!incindo o valor dos parâmetros de HSA ---
call atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,maxcont,0)
HM=> sai%HM(1:HMS)
ot => sai%ot


!_____________________________Inicializando HM ________________________

sai%HM(:)%fob = huge(0.d0) !Inicia a HM com valor enorme para a função objetivo

if(par%verb_t>0) write(*,'(A)')"HSA>> Formação da HM inicial: início"

do i = 1 , HMS !criando os pontos iniciais da HM
    
    do !obtendo um ponto aleatoriamente (verifica se é diferente dos anteriores)
        call rand_HSA_pto(par,HI)
        if(i==1) exit
        if( buscar_HSA_pto(sai%nAF, sai%hAF, HI) == 0 ) exit
    enddo
    
    !avaliando a função objetivo do ponto
    CALL fob(par%nvD, par%nvC, HI%vd, HI%vc, HI%fob )
    
    !armazenando o ponto com função objetivo avaliada
    sai%nAF = sai%nAF+1
    sai%hAF(sai%nAF) = HI
    call atualizar_HM(HI,HMS,HM,j) 
    
    !Definido ótimo na iteração inicial
    if(j==1) then !O novo ponto é o novo ótimo
        sai%ot_hAF(1) = sai%nAF
        sai%nAF_ot = sai%nAF
        sai%it_ot = 1
        ot = HI
    endif
enddo

if(par%verb_t>0) write(*,'(A)')"HSA>> Formação da HM inicial: fim"


!_____________________________bloco de iterações do algoritmo________________________

if(par%verb_t>0) write(*,'(A)')"HSA>> Iterações da otimização: início"

DO !loop das iterações no algoritmo de busaca pela harmonia 

    
    sai%it = sai%it + 1  !Incrementa o contador de iterações
    
    
    call atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,maxcont,cont)
    HM=> sai%HM(1:HMS)
   
    !improviso de uma nova harmonia ---- 
    
    !Laço sobre as variáveis discretas
    IF(par%nvD>0) THEN       
        DO i = 1 , par%nvD
            do j = 1 , HMS ;  HMvd(j) = HM(j)%vd(i) ; enddo !<-transfere os valores da vd(i) de todos os pontos da HM para HMvd(:)
            CALL improvisar_harmonia_vd(HMS, HMCRvd(i), PARvd(i), PAbw_vd(i), (/ par%Li_vD(i), par%Ls_vD(i) /), HMvd, HI%vd(i) )
        END DO 
    END IF
        
    !Laço sobre as variáveis contínuas
    IF(par%nvC>0) THEN       
        DO i = 1 , par%nvC
            do j = 1 , HMS ; HMvc(j) = HM(j)%vc(i) ; enddo !<-transfere os valores da vc(i) de todos os pontos da HM para HMvc(:)
           CALL improvisar_harmonia_vc(HMS, HMCRvc(i), PARvc(i), PAbw_vc(i), (/ par%Li_vC(i), par%Ls_vC(i) /), HMvc, HI%vc(i) )
        END DO 
            
    END IF
        

    !Avaliação da função objetivo ---
    
    !buscando por pontos avaliados anteriormente (apenas para problemas sem variáveis contínuas)
    i_hAF = 0
    if(par%nvC == 0) i_hAF = buscar_HSA_pto(sai%nAF, sai%hAF, HI)
         
    if(i_hAF/=0) then   !<-para de um ponto que já foi avaliado anteriormente
        HI%fob = sai%hAF(i_hAF)%fob
        sai%it_hAF(sai%it) = i_hAF
    
    else                !<-para de um ponto que não foi avaliado anteriormente
        CALL fob(par%nvD, par%nvC, HI%vd, HI%vc, HI%fob ) !<-avaliação da fob
        sai%nAF = sai%nAF+1
        sai%hAF(sai%nAF) = HI
        sai%it_hAF(sai%it) = sai%nAF

    endif
    
    
    !Atualiação da HM ---
    call atualizar_HM(HI,par%maxHMS,sai%HM,j)
    
    
    !Atualiza o valor ótimo ---
    novo_ot = j==1 !novo ótimo caso HI tenha sido adicionado no início da memória
    
    if(novo_ot) then !Atalização das variáveis em sai no caso de novo ótimo
        ot = HI
        sai%nAF_ot = sai%nAF
        sai%it_ot = sai%it
        sai%ot_hAF(sai%it) = sai%nAF
        
    elseif(sai%it>1) then !Manutenção de ot_hAF nessa iteração
        sai%ot_hAF(sai%it) = sai%ot_hAF(sai%it-1)
           
    endif
    
 
    !Imprimindo saída da iteração ---
1001 format("HSA>> cont="I6,"/",I6.6,2x,ES12.4E2)
    if(par%verb_t>0) write(*,1001) cont,maxcont,ot%fob
    
    !Verificar o encerramento das iterações do HSA  ---
    if(sai%it >= par%maxIT) exit
    if(sai%nAF >= par%maxAF) exit
    !if(cont>=maxcont) exit
    
    
    END DO !Fim do loop de iterações do HSA   

    call cpu_time(sai%tempo)
    sai%tempo = sai%tempo - tempo_inicio
 
    
    if(par%verb_t>0) then
        write(*,'(A)')"HSA>> Iterações da otimização: FIM"
        write(*,'(A,es12.4e2)')"HSA>> tempo da otimização: ",sai%tempo
        
        !    write(*,'(A)')"Resultado da otimização:"
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
!* PROPÓSITO:																		*
!* Gerar o valor para uma variável discreta na nova harmonia improvisada            *
!*                                            										*
!************************************************************************************
INTEGER,INTENT(IN):: HMS    !Tamnho da memoria de harmonia (Harmony Memory Size)
REAL(8),INTENT(IN):: HMCR   !Harmony Memory Considering Rate
REAL(8),INTENT(IN):: PAR    !Pitch adjusting rate
INTEGER,INTENT(IN):: PAbw   !Pitch adjusting distance bandwidth
INTEGER,INTENT(IN):: Lim(2) !Limites para o valor gerado (1-mínimo/2-máximo)
INTEGER,INTENT(IN):: HM(HMS)    !Memória de harmonia pra a variável em questão
INTEGER,INTENT(OUT):: x     !valor da variável na improvisação de um anova harmonia


!variáveis internas da subrotina
LOGICAL:: teste !auxiliar para decisões
INTEGER:: ix, ia    !auxiliar

!decidindo se o novo valor será gerado usando HM ou será gerado aleatoriamente
CALL lRAND(teste, HMCR)

IF(teste) THEN !gerar x usando HM

    !Escolhe aleatoriamente a posição em HM par usar na definição de x
    CALL iRAND(1,HMS,ix) 
    
    !Define o valor de x
    x = HM(ix)
 
    !decide se será realizado o ajuste (pitch adjustment)
    CALL lRAND(teste, PAR)
    
    IF(teste) THEN !realizar of pitch adjustment
        
        !gera o ajuste
        CALL iRAND(1,PAbw,ia)
        
        !decide a direção do ajuste na memória
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

END IF      !Fim do teste de utilização de HM na geração de x     

IA =0    
END SUBROUTINE improvisar_harmonia_vd




SUBROUTINE improvisar_harmonia_vc(HMS,HMCR,PAR,PAbw,Lim,HM,x)
!************************************************************************************
!* PROPÓSITO:																		*
!* Gerar o valor para uma variável discreta na nova harmonia improvisada            *
!*                                            										*
!************************************************************************************
INTEGER,INTENT(IN):: HMS    !Tamnho da memoria de harmonia (Harmony Memory Size)
REAL(8),INTENT(IN):: HMCR   !Harmony Memory Considering Rate
REAL(8),INTENT(IN):: PAR    !Pitch adjusting rate
REAL(8),INTENT(IN):: PAbw   !Pitch adjusting distance bandwidth
REAL(8),INTENT(IN):: Lim(2) !Limites para o valor gerado (1-mínimo/2-máximo)
REAL(8),INTENT(IN):: HM(HMS)    !Memória de harmonia pra a variável em questão
REAL(8),INTENT(OUT):: x     !valor da variável na improvisação de um anova harmonia


!variáveis internas da subrotina
LOGICAL:: teste !auxiliar para decisões
INTEGER:: ix    !auxiliar
REAL(8):: ra    !ajuste
REAL(8):: inf, sup  !valores inferior e superior que podem ser usados para PA

!decidindo se o novo valor será gerado usando HM ou será gerado aleatoriamente
CALL lRAND(teste, HMCR)

IF(teste) THEN !gerar x usando HM

    !Escolhe aleatoriamente a posição em HM par usar na definição de x
    CALL iRAND(1,HMS,ix) 
    
    !Define o valor de x
    x = HM(ix)
 
    !decide se será realizado o ajuste (pitch adjustment)
    CALL lRAND(teste, PAR)
    
    IF(teste) THEN !realizar of pitch adjustment
        
        !gera o ajuste
        CALL rRAND(0.d0,PAbw,ra)
        
        !decide a direção do ajuste na memória
        CALL lRAND(teste)
        
        IF(teste) THEN !ajuste para valor anterior
            x = max(Lim(1), x-ra )
        ELSE
            x = min(Lim(2), x+ra )
        END IF
              
    END IF          !fim do pitch adjustment

ELSE        !gerar x aleatoriamente

    CALL rRAND(lim(1),lim(2),x)

END IF      !Fim do teste de utilização de HM na geração de x     
    
END SUBROUTINE improvisar_harmonia_vc




subroutine atualizar_HM(pto,HMS,HM,i)
!************************************************************************************
!>Introduz um novo ponto "pto" na Harmony memory "HM" se o valor da função objetivo
!> do novo ponto for 'melhor' (menor valor) do que o pior componente da HM. 
!>A HM fica organizada em ordem crescente do valor da fob.
!>A melhor solução atual fica em HM(1) e a pior em HM(HMS)
!************************************************************************************
type(HSA_pto),intent(in)::    pto       !<Novo vetor solução
integer,intent(in):: HMS                !<Tamanho de HM
type(HSA_pto),intent(inout):: HM(:)     !<Harmony memory
integer,intent(out)::        i          !<Posição do novo vetor na HM (obs: i=0 se pto não for adicionado)

integer:: j, k
i=0 !Inicia com posição "0" indicando que a HM não foi atualizada.
if(pto%fob>HM(HMS)%fob) return !Sai caso pto seja pior que todos na HM
do j = 1 , HMS 
    if(pto%fob<HM(j)%fob) then !Encontrada a posição para inserir pto
        if(j/=HMS) then        !Deslocando pontos piores em uma posição
            do k = HMS , j+1 , -1 
                HM(k) = HM(k-1)
            enddo
        endif
        i=j                     !Posição de pto na HM
        HM(i)=pto               !Adicionando pto na HM
        exit                    !sai do laço quando encontrou a posição
    endif
enddo

end subroutine



!---------------------- rotinas relacionadas à variável tipo HSA_par ----------------



subroutine atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,imax,ii)
!************************************************************************************
!>Retorna o valor das variáveis com base no valor do contador atual
!************************************************************************************
type(HSA_par),intent(in):: par              !<Parâmetros informados inicialmente
integer,intent(out):: HMS                   !Tamanho da memória usada (HMS)
real(8),intent(out):: HMCRvd(:), PARvd(:)   !<HMCR e PAR para variáveis discretas
integer,intent(out):: PAbw_vd(:)            !bw para variáveis discretas
real(8),intent(out):: HMCRvc(:), PARvc(:)   !<HMCR e PAR par variáveis contínuas
real(8),intent(out):: PAbw_vc(:)            !<bw para variáveis contínuas
integer,intent(in):: imax                   !Valor limite do contador
integer,intent(in):: ii                     !<Contador

integer:: jj, tt
real(8):: rr, mm        !valor real de ii e imax, respectivamente
character(30) err_str !string com mensagem de erro

!variáveis para controle da impressão dos parâmetos em arquivo
logical:: arq_aberto !informação se o arquivo está aberto
character(50):: arq_action, frmt1

!----
300 format("ERR> @atualizar_parametros_HSA_par: em ii=",i9,A30)

    
err_str = ' '
rr = real(ii,8) !converte ii em real para ser usado na função fnc
mm = real(imax,8)

!---

HMS = nint( fnc(par%HMS_tf,par%HMS_c,mm,rr) )
if(HMS>par%maxHMS) then 
    write(err_str,"(' HMS = ',i3,' > maxHMS = ',i3)")HMS,par%maxHMS
    goto 301
endif

!---

if(par%nvD>0) then !parâmetros para as variáveis discretas
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

if(par%nvC>0) then !parâmetros para as variáveis contínuas
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

!!impressão do valor dos parâmetros em arquivo
call imprimir_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,ii,par%uni_imp_par)

!---

return !retorna sem problemas

301 write(*,300)ii,err_str !informa sobre o problema encontrado

end subroutine atualizar_parametros_HSA_par



subroutine imprimir_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,ii,uni)
!************************************************************************************
!>Imprime o valor das variáveis em um arquivo
!************************************************************************************
type(HSA_par),intent(in):: par              !<Parâmetros informados inicialmente
integer,intent(in):: HMS                   !Tamanho da memória usada (HMS)
real(8),intent(in):: HMCRvd(:), PARvd(:)   !<HMCR e PAR para variáveis discretas
integer,intent(in):: PAbw_vd(:)            !bw para variáveis discretas
real(8),intent(in):: HMCRvc(:), PARvc(:)   !<HMCR e PAR par variáveis contínuas
real(8),intent(in):: PAbw_vc(:)            !<bw para variáveis contínuas
integer,intent(in):: ii                     !<Contador
integer,intent(in):: uni                    !Unidade para a impressão dos parâmetros


integer:: jj, tt
character(30) err_str !string com mensagem de erro

!variáveis para controle da impressão dos parâmetos em arquivo
logical:: arq_aberto !informação se o arquivo está aberto
character(50):: arq_action, frmt1

!----
300 format("ERR> @imprimir_parametros_HSA_par: em ii=",i9,A30)

err_str = ' '

!impressão do valor dos parâmetros em arquivo
if(uni/=0) then

    if(ii==0) then !ações no início da otimização
        
        !verificando a conformidade do arquivo associado à unidade uni
        inquire(unit=uni, opened= arq_aberto,action=arq_action)
        if(.not. arq_aberto) then
            write(err_str,'(A,i5)') "Arquivo não aberto, unidade=" , uni
            goto 301
        else if(arq_action=='READ') then
            write(err_str,'(A,i5)') "Arquivo somente para leitura, unidade=" , uni
            goto 301
        endif        
        
        !cabeçalho com lista de variáveis (formato tecplot)  
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
    
    !imprimindo parâmetros
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
!>	Verifica se haverá problemas com os parâmetros gerados durante a otimização
!************************************************************************************
type(HSA_par),intent(in):: par              !<Parâmetros informados inicialmente
integer,intent(in):: imax                   !Valor limite do contador


integer:: HMS                       !Tamanho da memória usada (HMS)
real(8):: HMCRvd(maxNvD), PARvd(maxNvD)   !<HMCR e PAR para variáveis discretas
integer:: PAbw_vd(maxNvD)              !bw para variáveis discretas
real(8):: HMCRvc(maxNvC), PARvc(maxNvC)   !<HMCR e PAR par variáveis contínuas
real(8):: PAbw_vc(maxNvC)              !<bw para variáveis contínuas
integer:: ii                        !<Contador

character(1):: resp
integer:: uni , err
logical:: aberto

!Abrindo um arquivo para imprimir a verificação dos dados
uni=0
do
    uni = uni + 1
    inquire(unit=uni,opened=aberto)
    if(.not. aberto) open(unit=uni,file="verificar_parametos_HSA_par.plt",action="WRITE",iostat=err)
    if(err==0) exit
    if(err==500) then
        write(*,*)"ERR>> @verificar_parametos_HSA_par : não foi possível abrir um arquivo para imprimir &
        o arquivo de verificação"
        goto 401
    endif
    
enddo

if(par%verb_t>0) write(*,*)"@verificar_parametos_HSA_par : Verificando os parâmetros da otimização"

!Gerando e verificando os valores dos parâmetros que serão usados durante a otimização
do ii = 0 , imax
    call atualizar_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,imax,ii)
    call imprimir_parametros_HSA_par(par,HMS,HMCRvd,PARvd,PAbw_vd,HMCRvc,PARvc,PAbw_vc,ii,uni)
enddo

401 write(*,'(A)',advance='NO')"Prosseguir com a otimização? (S/N):"
read(*,*)resp
if(resp=='S' .or. resp=='s') then
    write(*,*) "Iniciando a otimização ... "
else
    write(*,*) "Abortando a otimização ... "
    stop
endif

close(uni)

end subroutine verificar_parametos_HSA_par



real(8) function fnc(tf,c,m,x)
!************************************************************************************
!>	Retorna o valor de uma das funções programadas 
!> O tipo de função é defindi por tf.
!> O valor é retornado para variável independente x, que deve estar entre 0 e m
!> m é o valor máximo de x
!************************************************************************************
integer,parameter:: nnc = ncoef     !<número de constantes
integer,parameter:: nntf = nctf     !<número de caracteres definindo a função
!Argumentos da função========================================
character(nntf),intent(in)::   tf   !<tipo de função
real(8),intent(in)::        c(nnc)  !<Coeficientes da função
real(8),intent(in)::        m       !<Valor máximo de x
real(8),intent(in)::        x       !<Variável independente

selectcase(tf) !selecionando o tipo de função
case('cst') !constante
    fnc = c(1)
case('lin') !linear ; fnc(0)=C1 e fnc(m)=C2
    fnc = c(1) + (c(2)-c(1))*x/m
case('pvf') !parábola com vértice no final; fnc(0)=c1 ; fnc(m)=c2 ; df(m)/dx = 0
    fnc = c(1) + (c(2)-c(1))* (x/m) * (2 - x/m) 
case('pv0') !parábola com vértice no início; fnc(0)=c1 ; fnc(m)=c2 ; df(0)/dx = 0
    fnc = c(1) + ( c(2)-c(1))*(x/m)**2
case('pim') !parábola com início e máximo/mínimo definido; fnc(0)=c1; fnc(c2)=c3; df(c2)/dx=0
    fnc = c(1) + (c(3)-c(1))* (x/c(2)) * (2 - x/c(2))
case('pg1') !polinômio do primeiro grau (linear)
    fnc = c(1) + c(2)*x
case('pg2') !polinômio do segundo grau (parábola)
    fnc = c(1) + c(2)*x + c(3)*x**2
case('pg3') !polinômio do terceiro grau (cúbica)
    fnc = c(1) + c(2)*x + c(3)*x**2 + c(4)*x**3
case('exp') !exponencial
    fnc = c(1) * exp( c(2) * x ) + c(3)
case default    
    fnc = huge(fnc)
    write(*,1001) tf
1001 format("ERR>> @fnc: tipo de funcao nao encontrada = ",A3)     
endselect

end function fnc



!---------------------- rotinas relacionadas à variável tipo HSA_pto ----------------

subroutine iniciar_HSA_pto(par,pto)
!************************************************************************************
!>Aloca as variáveis vd e vc de pto com tamanhos definidos em par%nvd e par%nvc
!>, além de zerar dotas as sub-variáveis de pto
!************************************************************************************
type(HSA_par),intent(in):: par  !<Parâmetros
type(HSA_pto),intent(out):: pto !<ponto no espaço de resposta

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
!> Desaloca as variáveis vd e vc e zera o valor de fob
!************************************************************************************
type(HSA_pto),intent(inout):: pto !<ponto no espaço de resposta
if(allocated(pto%vd)) deallocate(pto%vd)
if(allocated(pto%vc)) deallocate(pto%vc)
pto%fob = 0.d0
end subroutine terminar_HSA_pto



    
subroutine rand_HSA_pto(par,pto)
!************************************************************************************
!>Gera valores aleatórios para as variáves vd e vc, dentro dos limites estabelecidos
!> na variável par%Li_** e par%Li_** correspondente.
!************************************************************************************
type(HSA_par),intent(in):: par      !<Parâmetros
type(HSA_pto),intent(inout):: pto   !<ponto no espaço de resposta

if(par%nvD>0) call iRAND( par%Li_vD , par%Ls_vD, pto%vd) 
if(par%nvC>0) call rRAND( par%Li_vC , par%Ls_vC, pto%vc) 

end subroutine rand_HSA_pto    
    


function comparar_HSA_pto(pr,pc) result(comp)
!************************************************************************************
!> Compara dois HSA_pto, pc e pr, apresentando resutados em comp(i),
!> i=1 p/ var. discretas e i=2 para variáveis contínuas            
!*																					
!>	comp(i)= 0 - indica igualdade													
!>	comp(i)< 0 - indica que pc tem a variáveil de posição abs(comp(i)) menor que em pr		
!>	comp(i)> 0 - indica que pc tem a variáveil de posição comp(i) maior que em pr		
!>  comp(i)= huge - indica que pr e pc não tem o mesmo número de variáveis desse tipo
!> 	
!************************************************************************************
!parâmetros da subrotina
type(HSA_pto),intent(in):: pr   !<ponto de referência
type(HSA_pto),intent(in):: pc   !<ponto comparado com a referência
integer,dimension(2):: comp     !<vetor indicando a posição em que existe diferença

integer:: i, n!variáveis internas

comp = 0

!--- variáveis discretas
!Número de variáveis discertas 
n = size(pr%vd)
if(n/=size(pc%vd)) then !No caso de ter 
    comp(1) = huge(i)
else if(n>0) then      !No caso de haver variáveis discretas
    do i = 1 , n
        if(pc%vd(i)==pr%vd(i)) cycle
        comp(1) = sign(i, pc%vd(i)-pr%vd(i) )
        exit
    enddo
endif


!--- variáveis contínuas
!Número de variáveis contínuas em cada ponto
n = size(pr%vc)
if(n/=size(pc%vc)) then !No caso de ter 
    comp(2) = huge(i)
else if(n>0) then      !No caso de haver variáveis discretas
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
!< pos = 0 se pto não for encontrado em hAF 
!< pos igual à posição de pto em hAF							
!************************************************************************************
!parâmetros da subrotina
integer,intent(in):: nAF            !<número de pontos armazenados em hAF
type(HSA_pto),intent(in):: hAF(:)   !<pontos analisados anteriormente
type(HSA_pto),intent(in):: pto      !<ponto buscado

integer:: pos !<posição de pto em hAF
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
!>Impressão formatada da variável HSA_pto
!************************************************************************************
type(HSA_pto), intent(in) :: dtv
integer, intent(in) :: unit
character(len=*), intent(in) :: iotype
integer, intent(in) :: v_list(:)
integer, intent(inout) :: iostat
integer, intent(inout) :: iomsg

integer:: sl !número de informações passadas em v_list
integer:: nvd, lvd !número de variáveis discretas e núm. de dígitos para impressão
integer::nvc !Número de variáveis contínuas
integer:: f  !Código para o formato de impressão
character(50):: frmt !formato para impressão
integer:: ii    !variáveis auxiliares
 
sl = size(v_list)
lvd = 4 !Adotando variáveis inteiras com até 4 dígitos

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
    
    


!---------------------- rotinas relacionadas à variável tipo HSA_sai ----------------

    
subroutine iniciar_HSA_sai(par,sai)
!************************************************************************************
!* PROPÓSITO:																		*
!*																					*
!*																					*
!* PROGRAMADOR			: Felipe S. Almeida											*
!* Última alteração em	: agosto/2014												*
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
!>Desaloca as variáveis alocaveis de HSA_sai e zera as demais variáveis
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
!!>Imprime as informações sobre o processo de otimização em um arquivo
!!************************************************************************************
!type(HSA_par),intent(in):: par  !Parâmetros da otimização
!type(HSA_sai),intent(in)::  sai !Dados resultantes da otimização
!integer,intent(in):: uni        !Unidade para impressão dos dados contidos em sai
!
!!variáveis internas
!character(len=15):: rsp
!logical:: aberto
!integer:: ii, lvd
!character(300):: txt
!
!!----Verificando o arquivo atrelato à unidade uni
!inquire(unit=uni,opened=aberto,write=rsp)
!if(.not.aberto .or. rsp=="NO") then
!    write(*,*)"ERR>>@imprimir_HSA_sai: não é possível imprimir no arquivo"
!endif
!
!----
!ii = max( maxval( abs(par%Li_vd) ) , maxval( abs(par%Ls_vd)) )
!lvd = ceiling( log10( max( 9.d0 ,  real(ii , 8 ) ) ) ) ; pause
!
!end subroutine imprimir_HSA_sai




subroutine imprimir_HSA_sai_hAF_plt(par,sai,uni)
!************************************************************************************
!>Imprime as informações sobre o processo de otimização em um arquivo
!************************************************************************************
type(HSA_par),intent(in):: par  !Parâmetros da otimização
type(HSA_sai),intent(in)::  sai !Dados resultantes da otimização
integer,intent(in):: uni        !Unidade para impressão dos dados contidos em sai

!variáveis internas
integer:: jj, i, nit, j1,j2
character(300):: txt
character(100):: frmt

integer,allocatable:: aVd(:)
real(8),allocatable:: aVc(:)
real(8):: afob


write(uni,'(2/,A)') 'TITLE = "Historico HSA"'

!Linha com o nome das variáveis:
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