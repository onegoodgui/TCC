MODULE sislin_pbMEF

IMPLICIT NONE


TYPE MDGB_pbMEF
!****************************************************************************************
!* Variável tipo matriz banda geral (não simétrica) para uso nos programas de MEF      	*
!* Armazenagem dos dados conforme a biblioteca LAPACK pela nominação DGB				*
!* (D - real double precision) (GB - general band)  									*
!* OBS: A matriz é montada para representar uma matriz quadrada (nxn)      				*
!* A matriz é estruturalmente simétrica, ou seja:										*
!*	    Dada a matriz original A de dimensão nxn										*
!*		para cada coeficiente não nulo A(i,j) existe um coeficiente não nulo A(j,i)		*
!*		porém, em princípio A(i,j)/=A(j,i).                                      		*
!*	    Isso corresponde a uma matriz do tipo general band da LAPACK em que				*
!*	    o numero de diagonais superiores e inferiores (KU e KL, respecivamente)			*
!*		são iguais. Aqui KD=KI=ndsi														*
!*																						*
!****************************************************************************************
PRIVATE
    INTEGER:: LB                    !largura da banda
    INTEGER:: ntgl                  !Número total de graus de liberdade
    INTEGER:: ngln                  !Número de graus de liberdade por nó
    INTEGER:: nnos                  !Número de nós
    INTEGER:: ndsi                  !N. diagonais sup. ou inf. (iguais p/ matiz quadrad)
    INTEGER:: LD                    !variável LDAB da biblioteca LAPACK
    REAL(8),ALLOCATABLE:: M(:,:)    !Matriz com armazenamento banda geral
    INTEGER,ALLOCATABLE:: ipiv(:)   !índices p/ permutação de pivos (DGBSV-LAPACK)
END TYPE


TYPE MDSB_pbMEF
!****************************************************************************************
!* Variável tipo matriz banda simétrica para uso nos programas de MEF      	            *
!* Armazenagem dos dados conforme a biblioteca LAPACK pela nominação DSB				*
!* (D - real double precision) (SB - symmetric band)  									*
!* OBS: A matriz é montada armazenando a parte triangular superior da mantriz original  *
!*																						*
!****************************************************************************************
PRIVATE
    INTEGER:: LSB                   !largura da banda
    INTEGER:: ntgl                  !Número total de graus de liberdade
    INTEGER:: ngln                  !Número de graus de liberdade por nó
    INTEGER:: nnos                  !Número de nós
    INTEGER:: ndsi                  !N. diagonais sup. ou inf. (iguais p/ matiz quadrad)
    INTEGER:: LD                    !variável LDAB da biblioteca LAPACK
    REAL(8),ALLOCATABLE:: M(:,:)    !Matriz com armazenamento banda simétrica
END TYPE


!Definindo interface para algumas subrotinas definidas abaixo

!Inicializar (allocar) as matrizes banda (geral ou simétrica)
INTERFACE iniciar_MB_pbMEF
    MODULE PROCEDURE iniciar_MDGB_pbMEF, iniciar_MDSB_pbMEF
END INTERFACE iniciar_MB_pbMEF

!Zerar as matrizes banda (geral ou simétrica)
!Os dados das dimensões continuam salvos, apenas a matriz é zerada
INTERFACE zerar_MB_pbMEF
    MODULE PROCEDURE zerar_MDGB_pbMEF, zerar_MDSB_pbMEF
END INTERFACE zerar_MB_pbMEF

!Adiciona a matriz de um elemento à matrizes banda (geral ou simétrica)
INTERFACE adi_Me_MB_pbMEF
    MODULE PROCEDURE adi_Me_MDGB_pbMEF, adi_Me_MDSB_pbMEF
END INTERFACE adi_Me_MB_pbMEF

!Condição de contorno zero nas matrizes banda (geral ou simétrica)
!rotinas para aplicação do contorno em um único gdl (final _1)ou em um conjunto
!de gdl (final _n) indicado em um vetor "v_glc".
INTERFACE cont_zero_MB_pbMEF
    MODULE PROCEDURE    cont_zero_MDGB_pbMEF_1, cont_zero_MDGB_pbMEF_n, &
                        cont_zero_MDSB_pbMEF_1, cont_zero_MDSB_pbMEF_n
END INTERFACE cont_zero_MB_pbMEF

!Resolve o sistema de equações formado por uma matrizes banda (geral ou simétrica)
!e um array contendo NRHS vetores de termos independentes
INTERFACE resolver_SLE_MB_pbMEF
    MODULE PROCEDURE  resolver_SLE_MDGB_pbMEF, resolver_SLE_MDSB_pbMEF
END INTERFACE resolver_SLE_MB_pbMEF



CONTAINS

SUBROUTINE encontrar_MDC_pbMEF(nne,cone,MDC)
!****************************************************************************************
!* Calcuar a diferença máxima entre os nós conectados por um elemento e atualizar o 	*
!* valor MDC se for menor que a diferença encontrada no elemento						*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nne				    !Número de nós do elemento
INTEGER,INTENT(IN):: cone(nne)				!Conectividades do elemento
INTEGER,INTENT(INOUT):: MDC				    !Máxima diferença na conetividade da malha

!Variáveis internas
INTEGER:: cone_d(nne)   !Conetividades do elemento transladadas em um posição
INTEGER:: MDCe          !Máxima diferença na conetividade do elemento

cone_d(1:nne-1) = cone(2:nne)   ;   cone_d(nne) = cone(1) 

MDCe = MAXVAL( ABS( cone - cone_d ))

IF( MDCe > MDC ) MDC = MDCe

END SUBROUTINE encontrar_MDC_pbMEF



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!****************************************************************************************
! Subrotinas para matrizes banda geral (DGB)
!****************************************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SUBROUTINE iniciar_MDGB_pbMEF(nnos,ngln,MDC,MDGB,err)
!****************************************************************************************
!* Iniciar a variável do tipo MDGB_pbMEF                                             	*
!* Define o valor das variáveis auxiliares  											*
!* Faz a alocação e a inicialização das variáveis alocaveis						    	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nnos                   !Número de nós da malha
INTEGER,INTENT(IN):: ngln   				!Número de graus de liberdade por nó
INTEGER,INTENT(INOUT):: MDC				    !Máxima diferença na conetividade da malha
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   
INTEGER,INTENT(out):: err  				    !indicador de erro (err=0 ->sem erro)
err = 0
if(nnos<=0) err =1
if(ngln<=0) err =2
if(MDC <=0) err =3
if(err/=0) RETURN
MDGB%nnos = nnos 
MDGB%ngln = ngln
MDGB%ndsi = (MDC+1)*ngln-1  !Número de diagonais superiores ou inferiores
MDGB%LB = 2*MDGB%ndsi+1     !largura da banda
MDGB%LD = 3*MDGB%ndsi+1     !LDAB da LAPACK (conforme requisito da sub. DGBSV ) <<<<<<ver
MDGB%ntgl = nnos*ngln       !número total de graus de liberdade do sistema

ALLOCATE(MDGB%M(MDGB%LD,MDGB%ntgl),stat=err)    !allocando a matriz banda geral
if(err/=0) err = 4
MDGB%M = 0.d0

ALLOCATE(MDGB%ipiv(MDGB%ntgl),stat=err)    !allocando o vetor ipiv (DGBSV-LAPACK)
if(err/=0) err = 5
MDGB%ipiv = 0

END SUBROUTINE iniciar_MDGB_pbMEF


SUBROUTINE terminar_MDGB_pbMEF(MDGB,err)
!****************************************************************************************
!* Terminar a variável do tipo MDGB_pbMEF                                             	*
!* Zera o valor das variáveis auxiliares								    			*
!* Faz a desalocação das variáveis alocaveis            						    	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   
INTEGER,INTENT(out):: err  				    !indicador de erro (err=0 ->sem erro)
MDGB%nnos = 0 
MDGB%ngln = 0
MDGB%ndsi = 0  
MDGB%LB = 0    
MDGB%LD = 0    
MDGB%ntgl = 0  
DEALLOCATE(MDGB%M,stat=err)    !allocando a matriz banda geral
if(err/=0) err = 1
DEALLOCATE(MDGB%ipiv,stat=err)       !allocando o vetor ipiv (DGBSV-LAPACK)
if(err/=0) err = 2
END SUBROUTINE terminar_MDGB_pbMEF


SUBROUTINE zerar_MDGB_pbMEF(MDGB)
!****************************************************************************************
!* zera os coeficientes da matriz banda geral                                          	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   
MDGB%M = 0.d0
MDGB%ipiv = 0
END SUBROUTINE zerar_MDGB_pbMEF


SUBROUTINE somar_MDGB_pbMEF(C,A,B)
!****************************************************************************************
!* soma de duas matrizes banda geral                                                  	*
!*																						*
!****************************************************************************************
TYPE(MDGB_pbMEF),INTENT(IN)::  A        !matriz banda geral                   
TYPE(MDGB_pbMEF),INTENT(IN)::  B        !matriz banda geral                   
TYPE(MDGB_pbMEF),intent(inout)::  C     !matriz banda geral                   

INTEGER:: ndsi
LOGICAL:: err_dimen=.false.

!Conferindo erro nas dimensões das matrizes
IF(A%ntgl/=C%ntgl .or. B%ntgl /= C%ntgl) err_dimen = .true.
IF(A%LD/=C%LD .or. B%LD /= C%LD) err_dimen = .true.
IF(err_dimen) then
    write(*,*)"ERRO>>>somar_MDGB_pbMEF"
    read(*,'()')
    return
END IF

ndsi = C%ndsi
C%M(ndsi:,:) = A%M(ndsi:,:) + B%M(ndsi:,:)

END SUBROUTINE somar_MDGB_pbMEF


FUNCTION linha_GB(i,j,LB) result(iGB)
!****************************************************************************************
!* Determina a linha (iGB) na matriz armazenada em Banda Geral, correspondente a uma  	*
!* posição (i,j) em uma matriz  quadrada                 								*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: i				    !linha do coeficiente na matriz global padrao
INTEGER,INTENT(IN):: j				    !coluna do coeficiente na matriz global padrao
INTEGER,INTENT(IN):: LB				    !largura da banda
!resultado
INTEGER:: iGB				            !linha do coeficiente na matriz banda geral

iGB = LB+i-j
END FUNCTION


SUBROUTINE adi_Me_MDGB_pbMEF(nne,cone,Me,MDGB,err)
!****************************************************************************************
!* Adiciona a matriz de um elemento à matriz banda geral                               	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nne				    !Número de nós do elemento
INTEGER,INTENT(IN):: cone(nne)				!Conectividades do elemento
REAL(8),INTENT(IN):: Me(:,:)                !Matriz do elemento finito
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   
INTEGER,INTENT(OUT),OPTIONAL:: err          !indicador de erro

!variáveis internas
integer:: n, i, j, ii, jj, kk
integer:: LB                    !largura da banda
integer:: ngln                  !Múmero de gdl por nó
integer:: ndsi                  !Número de diagonais superiores ou inferiores
integer:: ngle                  !número de graus de liberdade do elemento
integer,allocatable:: gle(:)    !graus de liberdade globais do elemento       

LB = MDGB%LB
ngln = MDGB%ngln
ndsi = MDGB%ndsi      
ngle = nne*ngln       
allocate(gle(ngle))       

!calculando variáveis auxiliares para o apontamento dos componentes
!da matriz do elemento na matriz banda geral
ii = 0
do n = 1 , nne !laço sobre os nós do elemento
  kk = (cone(n)-1)*ngln 
  do j = 1 , ngln !laço sobre os gdl de cada nó
        ii = ii + 1
        gle(ii) = kk+j 
  end do
end do

IF(PRESENT(err)) THEN !verificar erros nos dados de entrada
    err = 0
    IF(MAXVAL(cone)>MDGB%nnos .or. MINVAL(cone)<=0 ) err = 1
    IF( (maxval(gle)-minval(gle)) > ndsi ) err = 2
    IF(.not. allocated(MDGB%M)) err = 3
    IF(err/=0) RETURN 
END IF

!somando os componentes da matriz do elemento na matriz global 
do j = 1 , ngle     !contador de colunas da matriz Me
    jj = gle(j)     !apontador da coluna em MDGB
    kk = LB - jj    !apontador da linha em MDGB  
    do i = 1 , ngle !contador de linhas da matriz Me
        ii = kk + gle(i)   !apontador da linha em MDGB 
        MDGB%M(ii,jj) = MDGB%M(ii,jj) + Me(i,j) 
    end do
end do
    
END SUBROUTINE adi_Me_MDGB_pbMEF



SUBROUTINE cont_zero_MDGB_pbMEF_1(glc,MDGB)
!****************************************************************************************
!* Aplica condição de contorno zero para uma variável (gl) na matriz banda geral       	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: glc				    !grau de liberdade com condição de contorno
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   

!variáveis internas
integer:: i, j, ii
integer:: ndsi                  !Número de diagonais superiores ou inferiores
integer:: li, ls                !limites inferior e superior p/ contador da coluna

!ndsi = MDGB%ndsi      
ndsi = MDGB%ndsi      

!zerando a coluna
MDGB%M(ndsi:,glc) = 0.d0        

!zerando a linha
li = max(1 , (glc - ndsi)) 
ls = min( MDGB%ntgl, (glc + ndsi) ) 
ii = MDGB%LB + glc

DO j = li , ls
    i = ii - j
    MDGB%M(i,j) = 0.d0
end do

!colocando valor 1 na diagonal principal
MDGB%M(MDGB%LB,glc) = 1.d0  
     
END SUBROUTINE cont_zero_MDGB_pbMEF_1


SUBROUTINE cont_zero_MDGB_pbMEF_n(nglc,v_glc,MDGB)
!****************************************************************************************
!* Aplica condição de contorno zero para uma variável (gl) na matriz banda geral       	*
!* Realiza a tarefa da mesma forma que na subrotina cont_zero_MDGB_pbMEF_1				*
!* Carrega grupo de glc (no vetor v_glc) p/ aplicar as condições de contorno na matriz. *
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nglc				    !Número de glc
INTEGER,INTENT(IN):: v_glc(nglc)			!grau de liberdade com condição de contorno
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   

!variáveis internas
integer:: n, glc            !<< adiconal em relação a cont_zero_MDGB_pbMEF_1 >>
integer:: i, j, ii
integer:: ndsi              !Número de diagonais superiores ou inferiores
integer:: li, ls            !limites inferior e superior p/ contador da coluna

!ndsi = MDGB%ndsi      
ndsi = MDGB%ndsi      

DO n = 1 , nglc !laço sobre os glc
    
    glc = v_glc(n)
    
    !início<<<<<trecho igual ao da cont_zero_MDGB_pbMEF_1 >>>>>
    
    !zerando a coluna
    MDGB%M(ndsi:,glc) = 0.d0        

    !zerando a linha
    li = max(1 , (glc - ndsi)) 
    ls = min( MDGB%ntgl, (glc + ndsi) ) 
    ii = MDGB%LB + glc

    DO j = li , ls
        i = ii - j
        MDGB%M(i,j) = 0.d0
    end do

    !colocando valor 1 na diagonal principal
    MDGB%M(MDGB%LB,glc) = 1.d0  
    
    !fim<<<<<trecho igual ao da cont_zero_MDGB_pbMEF_1 >>>>>

END DO
     
END SUBROUTINE cont_zero_MDGB_pbMEF_n


SUBROUTINE resolver_SLE_MDGB_pbMEF(MDGB,Y,err)
!****************************************************************************************
!* Resolve o sistema linear de equações A X = Y para X                                 	*
!* onde: A é a matriz MDGB, X é a matriz de solução e Y é a matriz de coeficientes		*
!* Y tem dimensões ntgl x NRHS                                               			*
!* X fica armazenado em Y na saida da subrotina		    								*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDGB_pbMEF),INTENT(INOUT)::  MDGB      !matriz banda geral                   
REAL(8),INTENT(INOUT):: Y(:,:)              !lado direito // resposta do sistema
INTEGER,INTENT(OUT),OPTIONAL:: err          !indicador de erro

!variáveis internas
INTEGER:: dim_Y(2)
EXTERNAL DGBSV
INTEGER:: err_acml

if(present(err)) err = 0 !<provisório

dim_Y = SHAPE(Y)

!Subrotina da LAPACK (chamada original)
!SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )

call DGBSV( MDGB%ntgl, MDGB%ndsi, MDGB%ndsi, dim_Y(2), MDGB%M, MDGB%LD, MDGB%ipiv, Y, dim_Y(1), err_acml )

END SUBROUTINE resolver_SLE_MDGB_pbMEF


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!****************************************************************************************
! Subrotinas para matrizes banda simétrica (DSB)
!****************************************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SUBROUTINE iniciar_MDSB_pbMEF(nnos,ngln,MDC,MDSB,err)
!****************************************************************************************
!* Iniciar a variável do tipo MDSB_pbMEF                                             	*
!* Define o valor das variáveis auxiliares  											*
!* Faz a alocação e a inicialização das variáveis alocaveis						    	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nnos                   !Número de nós da malha
INTEGER,INTENT(IN):: ngln   				!Número de graus de liberdade por nó
INTEGER,INTENT(INOUT):: MDC				    !Máxima diferença na conetividade da malha
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda simétrica                   
INTEGER,INTENT(out):: err  				    !indicador de erro (err=0 ->sem erro)
err = 0
if(nnos<=0) err =1
if(ngln<=0) err =2
if(MDC <=0) err =3
if(err/=0) RETURN
MDSB%nnos = nnos 
MDSB%ngln = ngln
MDSB%ndsi = (MDC+1)*ngln-1  !Número de diagonais superiores ou inferiores
MDSB%LSB =  MDSB%ndsi+1     !largura da semi-banda
MDSB%LD =   MDSB%LSB        !LDAB da LAPACK (conforme requisito da sub. DPBSV ) <<<<<<ver
MDSB%ntgl = nnos*ngln       !número total de graus de liberdade do sistema

ALLOCATE(MDSB%M(MDSB%LD,MDSB%ntgl),stat=err)    !allocando a matriz banda geral
if(err/=0) err = 4
MDSB%M = 0.d0

END SUBROUTINE iniciar_MDSB_pbMEF


SUBROUTINE terminar_MDSB_pbMEF(MDSB,err)
!****************************************************************************************
!* Terminar a variável do tipo MDSB_pbMEF                                             	*
!* Zera o valor das variáveis auxiliares								    			*
!* Faz a desalocação das variáveis alocaveis            						    	*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda geral                   
INTEGER,INTENT(out):: err  				    !indicador de erro (err=0 ->sem erro)
MDSB%nnos = 0 
MDSB%ngln = 0
MDSB%ndsi = 0  
MDSB%LSB = 0    
MDSB%LD = 0    
MDSB%ntgl = 0  
DEALLOCATE(MDSB%M,stat=err)    !allocando a matriz banda geral
if(err/=0) err = 1
END SUBROUTINE terminar_MDSB_pbMEF


SUBROUTINE somar_MDSB_pbMEF(C,A,B)
!****************************************************************************************
!* soma de duas matrizes banda simétrica                                              	*
!*																						*
!****************************************************************************************
TYPE(MDSB_pbMEF),INTENT(IN)::  A        !matriz banda simétrica                   
TYPE(MDSB_pbMEF),INTENT(IN)::  B        !matriz banda simétrica                   
TYPE(MDSB_pbMEF),INTENT(INOUT)::  C     !matriz banda simétrica                   

LOGICAL:: err_dimen=.false.

!Conferindo erro nas dimensões das matrizes
IF(A%ntgl/=C%ntgl .or. B%ntgl /= C%ntgl) err_dimen = .true.
IF(A%LD/=C%LD .or. B%LD /= C%LD) err_dimen = .true.
IF(err_dimen) then
    write(*,*)"ERRO>>>somar_MDSB_pbMEF"
    read(*,'()')
    return
END IF

C%M(:,:) = A%M(:,:) + B%M(:,:)

END SUBROUTINE somar_MDSB_pbMEF


SUBROUTINE zerar_MDSB_pbMEF(MDSB)
!****************************************************************************************
!* zera os coeficientes da matriz banda simétrica                                       *
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda simétrica                   
MDSB%M = 0.d0
END SUBROUTINE zerar_MDSB_pbMEF

FUNCTION linha_SB(i,j,LSB) result(iSB)
!****************************************************************************************
!* Determina a linha (iSB) na matriz armazenada em Banda Simétrica,  	                *
!* correspondente a uma posição (i,j) em uma matriz quadrada                 			*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: i				    !linha do coeficiente na matriz global padrao
INTEGER,INTENT(IN):: j				    !coluna do coeficiente na matriz global padrao
INTEGER,INTENT(IN):: LSB				!largura da semi-banda
!resultado
INTEGER:: iSB				            !linha do coeficiente na matriz banda geral

iSB = LSB+i-j
END FUNCTION


SUBROUTINE adi_Me_MDSB_pbMEF(nne,cone,Me,MDSB,err)
!****************************************************************************************
!* Adiciona a matriz de um elemento à matriz banda simétrica                            *
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nne				    !Número de nós do elemento
INTEGER,INTENT(IN):: cone(nne)				!Conectividades do elemento
REAL(8),INTENT(IN):: Me(:,:)                !Matriz do elemento finito
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda simétrica                   
INTEGER,INTENT(OUT),OPTIONAL:: err          !indicador de erro

!variáveis internas
integer:: n, i, j, ii, jj, kk
integer:: LSB                   !largura da semi-banda
integer:: ngln                  !Múmero de gdl por nó
integer:: ndsi                  !Número de diagonais superiores ou inferiores
integer:: ngle                  !número de graus de liberdade do elemento
integer,allocatable:: gle(:)    !graus de liberdade globais do elemento       

LSB = MDSB%LSB
ngln = MDSB%ngln
ndsi = MDSB%ndsi      
ngle = nne*ngln       
allocate(gle(ngle))       

!calculando variáveis auxiliares para o apontamento dos componentes
!da matriz do elemento na matriz banda geral
ii = 0
do n = 1 , nne !laço sobre os nós do elemento
  kk = (cone(n)-1)*ngln 
  do j = 1 , ngln !laço sobre os gdl de cada nó
        ii = ii + 1
        gle(ii) = kk+j 
  end do
end do

IF(PRESENT(err)) THEN !verificar erros nos dados de entrada
    err = 0
    IF(MAXVAL(cone)>MDSB%nnos .or. MINVAL(cone)<=0 ) err = 1
    IF( (maxval(gle)-minval(gle)) > ndsi ) err = 2
    IF(.not. allocated(MDSB%M)) err = 3
    IF(err/=0) RETURN 
END IF

!somando os componentes (triangular superior)
!da matriz do elemento na matriz global 
do j = 1 , ngle     !contador de colunas da matriz Me
    jj = gle(j)     !apontador da coluna em MDSB
    kk = LSB - jj   !apontador da linha em MDSB  
    do i = 1 , ngle !contador de linhas da matriz Me
        IF(gle(i)<=jj) THEN
            ii = kk + gle(i)   !apontador da linha em MDSB 
            MDSB%M(ii,jj) = MDSB%M(ii,jj) + Me(i,j)
        END IF 
    end do
end do
    
END SUBROUTINE adi_Me_MDSB_pbMEF


SUBROUTINE cont_zero_MDSB_pbMEF_1(glc,MDSB)
!****************************************************************************************
!* Aplica condição de contorno zero para uma variável (gl) na matriz banda simétrica    *
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: glc				    !grau de liberdade com condição de contorno
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda simétrica                   

!variáveis internas
integer:: i, j, ii
integer:: ndsi                  !Número de diagonais superiores ou inferiores
integer:: li, ls                !limites inferior e superior p/ contador da coluna

!ndsi = MDSB%ndsi      
ndsi = MDSB%ndsi      

!zerando a coluna
MDSB%M(:,glc) = 0.d0        

!zerando a linha
li = glc 
ls = min( MDSB%ntgl, (glc + ndsi) ) 
ii = MDSB%LSB + glc

DO j = li , ls
    i = ii - j
    MDSB%M(i,j) = 0.d0
end do

!colocando valor 1 na diagonal principal
MDSB%M(MDSB%LSB,glc) = 1.d0  
     
END SUBROUTINE cont_zero_MDSB_pbMEF_1


SUBROUTINE cont_zero_MDSB_pbMEF_n(nglc,v_glc,MDSB)
!****************************************************************************************
!* Aplica condição de contorno zero para uma variável (gl) na matriz banda simétrica    *
!* Realiza a tarefa da mesma forma que na subrotina cont_zero_MDSB_pbMEF_1				*
!* Carrega grupo de glc (no vetor v_glc) p/ aplicar as condições de contorno na matriz. *
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
INTEGER,INTENT(IN):: nglc				    !Número de glc
INTEGER,INTENT(IN):: v_glc(nglc)			!grau de liberdade com condição de contorno
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda geral                   

!variáveis internas
integer:: n, glc            !<< adiconal em relação a cont_zero_MDSB_pbMEF_1 >>
integer:: i, j, ii
integer:: ndsi              !Número de diagonais superiores ou inferiores
integer:: li, ls            !limites inferior e superior p/ contador da coluna

ndsi = MDSB%ndsi      

DO n = 1 , nglc !laço sobre os glc
    
    glc = v_glc(n)
    
    !início<<<<<trecho igual ao da cont_zero_MDSB_pbMEF_1 >>>>>
    
    !zerando a coluna
    MDSB%M(:,glc) = 0.d0        

    !zerando a linha
    li = glc 
    ls = min( MDSB%ntgl, (glc + ndsi) ) 
    ii = MDSB%LSB + glc

    DO j = li , ls
        i = ii - j
        MDSB%M(i,j) = 0.d0
    end do

    !colocando valor 1 na diagonal principal
    MDSB%M(MDSB%LSB,glc) = 1.d0  
    
    !fim<<<<<trecho igual ao da cont_zero_MDSB_pbMEF_1 >>>>>

END DO
     
END SUBROUTINE cont_zero_MDSB_pbMEF_n


SUBROUTINE resolver_SLE_MDSB_pbMEF(MDSB,Y,err)
!****************************************************************************************
!* Resolve o sistema linear de equações A X = Y para X                                 	*
!* onde: A é a matriz MDSB, X é a matriz de solução e Y é a matriz de coeficientes		*
!* Y tem dimensões ntgl x NRHS                                               			*
!* X fica armazenado em Y na saida da subrotina		    								*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda simétrica                   
REAL(8),INTENT(INOUT):: Y(:,:)              !lado direito // (saida=resposta do sistema)
INTEGER,INTENT(OUT),OPTIONAL:: err          !indicador de erro

!variáveis internas
INTEGER:: dim_Y(2)  !Dimensões do array Y
EXTERNAL DPBSV
INTEGER:: err_acml
INTEGER:: NRHS              !Múmero de sistemas de equações (number of right hand sides)

if(present(err)) err = 0 !  <provisório

dim_Y = SHAPE(Y)

NRHS = dim_Y(2)

!Subrotina da LAPACK (chamada original)
!SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

CALL DPBSV( 'U', MDSB%ntgl, MDSB%ndsi, NRHS, MDSB%M, MDSB%LD, Y, dim_Y(1), err_acml )

END SUBROUTINE resolver_SLE_MDSB_pbMEF


SUBROUTINE resolver_PAV_MDSB_pbMEF(MDSB,Y,err)
!****************************************************************************************
!* Resolve o problema de autovalores A X = lambda Y                                  	*
!* onde: A é a matriz MDSB, X é a matriz de solução e Y é a matriz de coeficientes		*
!* Y tem dimensões ntgl x NRHS                                               			*
!* X fica armazenado em Y na saida da subrotina		    								*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDSB_pbMEF),INTENT(INOUT)::  MDSB      !matriz banda simétrica                   
REAL(8),INTENT(INOUT):: Y(:,:)              !lado direito // (saida=resposta do sistema)
INTEGER,INTENT(OUT),OPTIONAL:: err          !indicador de erro

!variáveis internas
INTEGER:: dim_Y(2)  !Dimensões do array Y
EXTERNAL DPBSV
INTEGER:: err_acml
INTEGER:: NRHS              !Múmero de sistemas de equações (number of right hand sides)

if(present(err)) err = 0 !  <provisório

dim_Y = SHAPE(Y)

NRHS = dim_Y(2)

!Subrotina da LAPACK (chamada original)
!SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )

CALL DPBSV( 'U', MDSB%ntgl, MDSB%ndsi, NRHS, MDSB%M, MDSB%LD, Y, dim_Y(1), err_acml )

END SUBROUTINE resolver_PAV_MDSB_pbMEF



SUBROUTINE resolver_PAVG_MDSB_pbMEF(A,B,nRES,nAV,aVal,aVet,err)
!****************************************************************************************
!* Resolve o problema de autovalores generalizado Ax = lambda B x                       *
!* São calculados os primeiros nAV autovalores e autovetores                    		*
!*																						*
!****************************************************************************************
!Parâmetros da subrotina
TYPE(MDSB_pbMEF),INTENT(INOUT)::  A         !matriz banda simétrica A                   
TYPE(MDSB_pbMEF),INTENT(INOUT)::  B         !matriz banda simétrica B                  
integer,intent(in):: nRES                   !Número de gdls restritos
integer,intent(in):: nAV                    !Número de autovalores a serem calculados
real(8),intent(out):: aVal(:)               !Autovalores calculados
real(8),intent(out):: aVet(:,:)             !Autovetores calculados
INTEGER,INTENT(OUT),OPTIONAL:: err          !indicador de erro


!Variáveis internas
!Mesmos nomes usados na LAPACK,
!a menos de LD_ que substitui LDQ, LDZ e N em muitos casos
integer:: LD_ , M, info, i, j, jj, NSPLIT
integer,allocatable:: ifail(:), IWORK(:)
real(8):: ABSTOL
real(8),allocatable:: Q(:,:), W(:), Z(:,:), WORK(:) 


!Teste das dimensões das variáveis de entrada

LD_ = A%ntgl
allocate( Q(LD_,LD_) ) ; Q =0.d0
allocate( W(LD_) )     ; W = 0.d0
allocate( Z(LD_,LD_) ) ; Z =0.d0
allocate( WORK(7*LD_)) ; Work = 0.d0
allocate( iWORK(5*LD_) ) ; iwork = 0
allocate(ifail(nAV))     ; ifail = 0

!SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Q, LDQ,
                  ! VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )

ABSTOL = 2.d0 * dlamch__()
!write(*,*) "abstol=", abstol

call DSBGVX('V', 'I', 'U', LD_, A%ndsi, B%ndsi, A%M, A%LD, B%M, B%LD, Q, LD_, &
!             0.d0, 0.d0, 1 , nAV, ABSTOL, M, W, Z, LD_, WORK, iWORK, ifail, info     )
             0.d0, 0.d0, nRES+1, nRES+nAV, ABSTOL, M, W, Z, LD_, WORK, iWORK, ifail, info     )


!mensagens de erro
if(info>0) then  
    if(info<=nAV) then
        do i = 1 , info
            write(*,'(2/,A,I9)')"ERRO>> resolver_PAVG_MDSB_pbMEF :: falha ao calcular o autovetor/autovalor num.", ifail(i)
        enddo
    else
        write(*,'(2/,A,I9)')"ERRO>> resolver_PAVG_MDSB_pbMEF :: Falha ao fatorar matriz B (nao e positivo definida)",info - LD_
    endif
elseif(info<0) then
    write(*,'(2/,A,I9)')"ERRO>> resolver_PAVG_MDSB_pbMEF :: argumento invalido na Bib LAPACK.", info
endif

if(present(err)) err=info

!verificando autovares e autovetores correspondentes a gdl restritos
jj = 0
do i = 1 , nAV
    if(W(i)==1.d0) then
        if(maxval(abs(z(1:LD_,i)))==1.d0) then
            jj = jj + 1
            !write(*,*)"resolver_PAVG_MDSB_pbMEF>>> autovetor num.", i," nulo"
        endif
    endif
enddo


!colocando os valores obtidos
aVal(1:nAV)         = W(1:nAV)
aVet(1:LD_,1:nAV)   = Z(1:LD_,1:nAV)


!tentando calcular os JJ primeiros autovalores
if(jj/=0) then

   !   CALL DSBGST( 'V', 'U', LD_, A%ndsi, B%ndsi, A%M, A%LD, B%M, B%LD, Q, LD_, &
   !               WORK, INFO )


    !CALL DSBTRD( 'U', 'U', LD_, A%ndsi, A%M, A%LD, WORK( 1:LD_ ), &
    !               WORK( LD_+1:2*LD_ ), Q, LD_, WORK( 2*LD_+1: ), INFO )
    
    !obs: refazer os componentes da matriz tridiagonal ( WORK( 1:LD_ ) e  WORK( LD_+1:2*LD_ ) )
    !      usando os componentes na matriz A .... ver DSBTRD

    work = 0.d0
    do i = 1 , LD_
        work(i) = A%M(A%ndsi+1,i)
    enddo
    do i = 1, LD_ - 1
        work(i+LD_) = A%M(A%ndsi,i+1)
    enddo

    !cálculo dos autovalores
    CALL DSTEBZ( 'I', 'B', LD_, 0.d0, 0.d0, 1, JJ , ABSTOL, &
                    WORK( 1 ), WORK( LD_+1 ), M, NSPLIT, W, &
                   IWORK( 1 ), IWORK( LD_+1 ), WORK( 2*LD_+1 ), &
                   IWORK( 2*LD_+1 ), INFO )

    !cálculo dos autovetores
    CALL DSTEIN( LD_, WORK( 1 ), WORK( LD_+1 ), M, W, &
                    IWORK( 1 ), IWORK( LD_+1 ), Z, LD_, &
                    WORK( 2*LD_+1 ), IWORK( 2*LD_+1 ), IFAIL, INFO ) 

         DO J = 1, M
            CALL DCOPY( LD_, Z( 1 , J ), 1, WORK( 1 ), 1 )
            CALL DGEMV( 'N', LD_, LD_, 1.d0, Q, LD_, WORK, 1, 0.d0, &
                       Z( 1 , J ), 1 )
        enddo

    !colocando os valores obtidos
    aVal(1:JJ)         = W(1:JJ)
    aVet(1:LD_,1:JJ)   = Z(1:LD_,1:JJ)

endif

deallocate(Q)
deallocate(W)
deallocate(Z)
deallocate( WORK  )
deallocate( iWORK )
deallocate( ifail )

END SUBROUTINE resolver_PAVG_MDSB_pbMEF


real(8) function dlamch__()
!****************************************************************************************
!* Retorna um valor pequeno para tolerância de erro                                     *
!* OBS:: substitui a função auxiliar dlamch() da lapack no caso de argumento 'S'        *
!* que representa o valor safe minimum													*
!*																						*
!****************************************************************************************
real(8):: rnd, eps, sfmin, small, rmach
real(8),parameter:: one = 1.0d+0 , zero = 0.0d+0

rnd = one

if( one.eq.rnd ) then
    eps = epsilon(zero) * 0.5
else
    eps = epsilon(zero)
end if

sfmin = tiny(zero)
small = one / huge(zero)
if( small.ge.sfmin ) then
    sfmin = small*( one+eps )
end if
dlamch__ = sfmin

end function dlamch__

END MODULE sislin_pbMEF

