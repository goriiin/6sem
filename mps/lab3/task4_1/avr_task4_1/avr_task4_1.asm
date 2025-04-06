.include "m8515def.inc" ; ���� ����������� ��� ATmega8515
.def temp = r16 ;��������� �������
.equ START = 0 ;0-� ��� ����� PA
.equ STOP = 1 ;1-� ��� ����� PA
;***������� ��� ��������� �������/��������
.equ PRESCALER_64 = ((1<<CS01)|(1<<CS00))
.equ TIMER_SETTINGS = ((1<<WGM01)|(1<<COM00)|PRESCALER_64)
.equ VALUE = 109 ; ������� 0,5 ���
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
	out PORTB,temp ;������ ������� ������
	clr temp ;�������������
	out DDRA,temp ; ����� PA �� ����
	ldi temp,0b00000011 ;��������� �������������
	out PORTA,temp ; ���������� ����� PA
	ldi temp,VALUE ;��������� ���������
	out OCR0,temp ; �������� �����
test_start:
	sbic PINA,START ;�������� ���������
	rjmp test_stop ; ������ START
	ldi temp,TIMER_SETTINGS
	out TCCR0,temp ; ��������� �������
wait_0:
	sbis PINA,START ;�������� ����������
	rjmp wait_0 ; ������
test_stop:
	sbic PINA,STOP ;�������� ���������
	rjmp test_start ; ������ stop
	clr temp
	out TCCR0,temp ;���������� �������
wait_1:
	sbis PINA,STOP ;�������� ����������
	rjmp wait_1 ; ������
	rjmp test_start
