   M S B A M L   `   `   `  ����   	WPF3dTube    WPF3dTube.MainWindow   \ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     -     � Title$  
MainWindow��$	��350��$	��525��.��� 
������� � RowDefinitions �� $
��Auto ��� $��* �.+�'� -    '�Name$ viewport��$��1�� '� Camera =� -   
camera�� � FarPlaneDistance$ 50Q� � NearPlaneDistance$ 0Q� � LookDirection$ 0,0,-105� � UpDirection$ 0,1,05� � Position$ 0,0,301�	 =� FieldOfView$	 45Q�x�g� -   	model��
 g� Content
 h� -   	group��.����  �� Color$ White��t� $ White�� t� 	Direction$ 0 -1 -15��� -   	spotLight��$ #666666�� �� 	Direction$ 0 -1 -15� �� InnerConeAngle$ 30Q� �� OuterConeAngle$ 60Q� !� Position$ 0 1 301� � Triggers :�  :� RoutedEvent$ 
MouseEnter��.[��� .���� -    g�Name$ myStoryBoardY��