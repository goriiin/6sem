.include "m8515def.inc" ; ���� ����������� ��� ATmega8515

.def temp = r16 ;��������� �������

.equ G = 0 ;0-� ��� ����� PA
.equ D = 1 ;1-� ��� ����� PA
.equ A = 2 ;1-� ��� ����� PA
.equ F = 3 ;1-� ��� ����� PA


;***������� ��� ��������� �������/��������
.equ PRESCALER_64 = ((1<<CS01)|(1<<CS00))

.equ TIMER_SETTINGS = ((1<<WGM01)|(1<<COM00)|PRESCALER_64)

.equ VALUE_G = 146 ; ������� 0,5 ���
.equ VALUE_D = 195 ; ������� 0,5 ���
.equ VALUE_A = 130 ; ������� 0,5 ���
.equ VALUE_F = 164 ; ������� 0,5 ���


.macro timer_on
	sbi  PORTB, 0 
	ldi temp,TIMER_SETTINGS 
    out TCCR0,temp 
.endmacro 

.macro timer_off
	clr temp 
	out TCCR0,temp
.endmacro 



;***������� �������� ����������

.org $000
	 rjmp INIT
;***������������� ��

INIT:
	ldi temp,high(RAMEND) ;���������
	out SPH,temp ; ��������� �����
	ldi temp,low(RAMEND) ; �� ���������
	out SPL,temp ; ������ ��
	ser temp ;������������� �������

	out DDRB,temp ; ����� PB �� �����
	cbi PORTB, 0
	;out PORTB,temp ;������ ������� ������

	clr temp ;�������������
	out DDRA,temp ; ����� PA �� ����
	ldi temp,0b00001111 ;��������� �������������
	out PORTA,temp ; ���������� ����� PA
	

START:
	timer_off
check_G:
	sbic PINA, G
	rjmp check_D
; ***** �������� ������� ���� ��� ���� G 196��  ���� 146
	ldi temp,VALUE_G ;��������� ���������
	out OCR0,temp ; �������� �����
	timer_on

	play_G:
		sbis PINA, G
		rjmp play_G
	
	timer_off


check_D:
	sbic PINA, D
	rjmp check_A
; ***** �������� ���� ��� ���� D 146,8��  -- 195
	ldi temp,VALUE_D ;��������� ���������
	out OCR0,temp ; �������� �����
	
	timer_on
	
	play_D:
		sbis PINA, D
		rjmp play_D
	timer_off
check_A:
	sbic PINA, A
	rjmp check_F
; ***** �������� ���� ��� ���� A 220��  -- 130
	ldi temp,VALUE_A ;��������� ���������
	out OCR0,temp ; �������� �����
	

	timer_on

	play_A:
		sbis PINA, A
		rjmp play_A

	timer_off
check_F:
	sbic PINA, F
	rjmp START
	; ***** �������� ���� ��� ���� A 174,6��  -- 164
	ldi temp,VALUE_F ;��������� ���������
	out OCR0,temp ; �������� �����

	timer_on

	play_F:
	    sbis PINA, F
	    rjmp play_F
	    
	    timer_off
	    rjmp START
