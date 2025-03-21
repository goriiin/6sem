.include "m8515def.inc" ; ����������� ATmega8515

.def temp = r16 ; ��������� �������

;*** ������� ��������
.org $000
	rjmp INIT ; ��������� ������
	rjmp PROGRAM_CALL

.org $007
	rjmp T0_OVF ;��������� ������� T0


;*** ������������� ��
INIT:
	ldi temp,low(RAMEND) ; ��������� ��������� �����
	out SPL,temp 
	ldi temp,high(RAMEND) 
	out SPH,temp 

	clr temp ;  PB �� ����
	out DDRB,temp ;
    

	; (1<<PB0)

	ldi temp,(1<<PB2) ; ��������� �������������� 
	out PORTB,temp ; ��������� ����� PB7

	ser temp ; ������������� ������� ����� PD
	out DDRD,temp ; �� �����
	out PORTD,temp ; ���������� �����������

	ldi temp,(1<<SE) ; ���������� ��������
	out MCUCR,temp ; � ����� Idle ����� ������ sleep

	;***��������� ������� T0 �� ����� �������� �������
	ldi temp,0x02 ; ���������� ����������
	out TIMSK,temp ; �� ������������ T0

	ldi temp,0x07 ; ������������ �������
	out TCCR0,temp ; �� �������������� �������� ����������

	sei ;���������� ���������� ����������
	ldi temp,0xFC ;$ FC=-4 ���
	out TCNT0,temp ; ������� 4� �������


LOOP:
	sleep ; ������� � ����� ����������� ����������������� 
	nop ; 
	rjmp LOOP

	;*** ��������� ���������� ��� ������������ ������� �0
T0_OVF:
	clr temp
	out PORTD,temp ; ��������� �����������

	rcall DELAY ; ��������

	ser temp
	out PORTD,temp ; ���������� �����������

	ldi temp,0xFC ; ������������
	out TCNT0,temp ; TCNT0

	reti


PROGRAM_CALL:
    ; �������� ��������� ������ � ���������
    sbic PORTB, PB0
    rjmp PROGRAM_CALL ; ���� ������ �� ������, ���������� ��������

    ; ��������� �������� ��� ���������� ��������
    rcall DELAY

    ; ����������� ��������� �������
    cbi PORTB, 0
    sbi PORTB, 0

;*** �������� ***
DELAY:
	ldi r19,6
	ldi r20,255
	ldi r21,255

dd:
	dec r21
	brne dd
	
	dec r20
	brne dd
	
	dec r19
	brne dd
	
	ret

;*** �������� ***
DELAY2:
	ldi r20,100
	ldi r21,255

dd:
	dec r21
	brne dd
	
	dec r20
	brne dd
	
	ret

